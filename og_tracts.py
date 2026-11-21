#!/usr/bin/env python3
"""Step 2c — attach tracts to orthogroups and pick the structure set.

For every 1:1 orthogroup, record per species whether the protein carries an
fLPS region (any class, and by class) and its longest poly-N run. Then select
the orthogroups for structural analysis:

  plasmodium_complete == 1 (one ortholog in every Plasmodium species), and
  at least one Plasmodium member carries a compositionally biased region.

Every member of a selected group is modelled, including members WITHOUT a
region: they are the reference that an insertion is mapped against.

Outputs
  og_tracts.tsv         one row per 1:1 orthogroup, per-species tract columns
  structure_set.tsv     species, gene, uniprot, og_id for every model to fetch
  og_tracts_summary.md  counts, including insertion presence/absence patterns
"""
import argparse
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polylib as P  # noqa: E402

CLASSES = ["N-rich", "GS-rich", "charged", "Q-rich", "other"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--orthogroups", required=True)
    ap.add_argument("--lcr-dir", required=True)
    ap.add_argument("--seqs", required=True, help="unmasked <key>_proteins.fasta")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    from call_tracts import load_proteins, max_run  # noqa: E402

    species = [r["key"] for r in P.read_tsv(a.species)]
    plas = [s for s in species if s.startswith("p")]

    lcr = {s: defaultdict(Counter) for s in species}
    lcr_res = {s: Counter() for s in species}
    for s in species:
        path = os.path.join(a.lcr_dir, f"{s}_lcr.tsv")
        if os.path.exists(path):
            for r in P.read_tsv(path):
                lcr[s][r["gene"]][r["lcr_class"]] += 1
                lcr_res[s][r["gene"]] += int(r["length"])
    seqs = {}
    for s in species:
        prots, _ = load_proteins(os.path.join(a.seqs, f"{s}_proteins.fasta"))
        seqs[s] = {g: seq for g, (_, seq) in prots.items()}

    rows, chosen = [], []
    pattern = Counter()
    class_in_pf_only = Counter()
    for og in P.read_tsv(a.orthogroups):
        if og["status"] != "1to1":
            continue
        row = {"og_id": og["og_id"], "n_species": og["n_species"],
               "plasmodium_complete": og["plasmodium_complete"]}
        has = {}
        for s in species:
            g = og[f"{s}_gene"]
            c = lcr[s].get(g, Counter()) if g else Counter()
            seq = seqs[s].get(g, "") if g else ""
            row[f"{s}_gene"] = g
            row[f"{s}_uniprot"] = og.get(f"{s}_uniprot", "")
            row[f"{s}_length"] = len(seq)
            row[f"{s}_lcr_n"] = sum(c.values())
            row[f"{s}_lcr_res"] = lcr_res[s].get(g, 0)
            row[f"{s}_lcr_classes"] = ",".join(f"{k}:{v}" for k, v in sorted(c.items()))
            row[f"{s}_max_n_run"] = max_run(seq, "N") if seq else ""
            has[s] = bool(c)
        row["plasmodium_lcr_pattern"] = "".join("1" if has.get(s) else "0" for s in plas)
        rows.append(row)
        if og["plasmodium_complete"] == "1":
            pattern[row["plasmodium_lcr_pattern"]] += 1
            if any(has[s] for s in plas):
                for s in species:
                    if row[f"{s}_gene"]:
                        chosen.append({"og_id": og["og_id"], "species": s,
                                       "gene": row[f"{s}_gene"],
                                       "uniprot": row[f"{s}_uniprot"],
                                       "has_lcr": int(has[s])})

    fields = ["og_id", "n_species", "plasmodium_complete", "plasmodium_lcr_pattern"]
    for s in species:
        fields += [f"{s}_{x}" for x in ("gene", "uniprot", "length", "lcr_n", "lcr_res",
                                        "lcr_classes", "max_n_run")]
    os.makedirs(a.out, exist_ok=True)
    P.write_tsv(os.path.join(a.out, "og_tracts.tsv"), rows, fields)
    P.write_tsv(os.path.join(a.out, "structure_set.tsv"), chosen,
                ["og_id", "species", "gene", "uniprot", "has_lcr"])

    complete = sum(pattern.values())
    sel_ogs = len({c["og_id"] for c in chosen})
    by_sp = Counter(c["species"] for c in chosen)
    mapped = Counter(c["species"] for c in chosen if c["uniprot"])
    lines = [
        "# Step 2c — tracts on orthogroups",
        "",
        f"1:1 orthogroups with one ortholog in every *Plasmodium* species: **{complete:,}**",
        f"of which at least one *Plasmodium* member carries an fLPS region: **{sel_ogs:,}**",
        "",
        "## Region presence across *Plasmodium* orthologs",
        "",
        f"Pattern order: {', '.join(plas)} (1 = protein carries a region).",
        "",
        "| pattern | orthogroups |",
        "|---|---:|",
    ]
    for pat, n in pattern.most_common():
        lines.append(f"| `{pat}` | {n:,} |")
    lines += [
        "",
        "## Structure set",
        "",
        "| species | models | with UniProt id |",
        "|---|---:|---:|",
    ]
    for s in species:
        lines.append(f"| {s} | {by_sp[s]:,} | {mapped[s]:,} |")
    lines.append(f"| **total** | **{sum(by_sp.values()):,}** | **{sum(mapped.values()):,}** |")
    lines += ["", f"At ~360 KB per AlphaFold PDB file, about "
                  f"{sum(mapped.values()) * 0.36 / 1024:.1f} GB to download.", ""]
    with open(os.path.join(a.out, "og_tracts_summary.md"), "w") as fh:
        fh.write("\n".join(lines))
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
