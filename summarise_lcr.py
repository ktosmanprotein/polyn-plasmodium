#!/usr/bin/env python3
"""Step 2a — compositionally biased regions with fLPS 2 (Harrison 2017, 2021).

Replaces the SEG-like caller of step 1 for region classes. fLPS tests each
window against a background composition with a binomial P-value and reports
which residues drive the bias, so the region class comes from the tool's own
signature rather than from ad hoc composition rules.

Background is each species' OWN proteome composition, so a region is called
when it is biased relative to that proteome. This matters across species with
very different AT content: against a generic background, AT-rich proteomes
look biased almost everywhere (tested: 94% of P. falciparum residues).

Input : <key>.fasta (all isoforms) and fLPS long-format output for it.
Output: <key>_lcr.tsv, one row per merged region on the primary isoform,
        lcr_summary.tsv across species.
"""
import argparse
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polylib as P  # noqa: E402

LCR_FIELDS = ["gene", "uniprot", "start", "end", "length", "lcr_class",
              "signatures", "min_p", "rel_position", "n_frac", "gs_frac", "ked_frac"]
SUMMARY_FIELDS = ["species", "proteins", "proteins_with_lcr_pct", "lcr_residue_pct",
                  "regions", "N-rich", "GS-rich", "charged", "Q-rich", "other",
                  "proteins_with_N-rich_pct", "proteins_with_GS-rich_pct"]


def signature_class(sig):
    """Class of one fLPS signature such as '{N}', '{N,K}', '{S,G}'.

    The first residue is the dominant one (fLPS lists them in order of
    precedence)."""
    res = [r for r in sig.strip("{} ").split(",") if r]
    if not res:
        return "other"
    s = set(res)
    if res[0] == "N":
        return "N-rich"
    if res[0] == "Q":
        return "Q-rich"
    if s <= set("GS"):
        return "GS-rich"
    if s <= set("KEDR"):
        return "charged"
    return "other"


CLASS_PRIORITY = ["N-rich", "Q-rich", "GS-rich", "charged", "other"]


def parse_flps(path):
    """Yield (seq_id, start, end, signature, p) for SINGLE/MULTIPLE biases;
    coordinates converted to 0-based half-open."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 8 or p[1] not in ("SINGLE", "MULTIPLE"):
                continue
            try:
                yield p[0], int(p[3]) - 1, int(p[4]), p[7], float(p[6])
            except ValueError:
                continue


def merge_regions(hits):
    """Merge overlapping (start, end, sig, p) hits into regions; each region
    keeps all signatures and the smallest P-value."""
    out = []
    for s, e, sig, p in sorted(hits):
        if out and s < out[-1]["end"]:
            r = out[-1]
            r["end"] = max(r["end"], e)
            r["sigs"].append(sig)
            r["p"] = min(r["p"], p)
        else:
            out.append({"start": s, "end": e, "sigs": [sig], "p": p})
    return out


def region_class(sigs):
    """Class of a merged region: the highest-priority class among its
    SINGLE signatures if any, else among all signatures."""
    classes = {signature_class(s) for s in sigs}
    for c in CLASS_PRIORITY:
        if c in classes:
            return c
    return "other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--seqs", required=True)
    ap.add_argument("--flps-dir", required=True)
    ap.add_argument("--apico-proteomes", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from call_tracts import uniprot_map, load_proteins  # noqa: E402

    os.makedirs(a.out, exist_ok=True)
    summary = []
    for sp in P.read_tsv(a.species):
        key = sp["key"]
        flps = os.path.join(a.flps_dir, f"{key}_flps.tsv")
        prot = os.path.join(a.seqs, f"{key}_proteins.fasta")
        if not (os.path.exists(flps) and os.path.exists(prot)):
            print(f"[skip] {key}: missing {flps if not os.path.exists(flps) else prot}")
            continue
        proteins, _ = load_proteins(prot)
        id2gene = {m["id"]: g for g, (m, _) in proteins.items()}
        to_uni = uniprot_map(a.apico_proteomes, key)

        hits = defaultdict(list)
        for sid, s, e, sig, p in parse_flps(flps):
            g = id2gene.get(sid)
            if g:
                hits[g].append((s, e, sig, p))

        rows, cls_n, cls_prot = [], Counter(), defaultdict(set)
        lcr_res = 0
        for g, h in sorted(hits.items()):
            seq = proteins[g][1]
            L = len(seq)
            for r in merge_regions(h):
                seg = seq[r["start"]:r["end"]]
                c = region_class(r["sigs"])
                cls_n[c] += 1
                cls_prot[c].add(g)
                lcr_res += r["end"] - r["start"]
                n = max(len(seg), 1)
                rows.append({
                    "gene": g, "uniprot": to_uni.get(g, ""),
                    "start": r["start"] + 1, "end": r["end"],
                    "length": r["end"] - r["start"], "lcr_class": c,
                    "signatures": ";".join(dict.fromkeys(r["sigs"])),
                    "min_p": f"{r['p']:.2e}",
                    "rel_position": f"{((r['start'] + r['end']) / 2) / L:.3f}",
                    "n_frac": f"{seg.count('N') / n:.3f}",
                    "gs_frac": f"{(seg.count('G') + seg.count('S')) / n:.3f}",
                    "ked_frac": f"{sum(seg.count(x) for x in 'KED') / n:.3f}",
                })
        P.write_tsv(os.path.join(a.out, f"{key}_lcr.tsv"), rows, LCR_FIELDS)
        n_prot = len(proteins)
        total = sum(len(s) for _, s in proteins.values())
        summary.append({
            "species": key, "proteins": n_prot,
            "proteins_with_lcr_pct": f"{100 * len(hits) / n_prot:.1f}",
            "lcr_residue_pct": f"{100 * lcr_res / total:.1f}",
            "regions": len(rows),
            **{c: cls_n[c] for c in CLASS_PRIORITY},
            "proteins_with_N-rich_pct": f"{100 * len(cls_prot['N-rich']) / n_prot:.1f}",
            "proteins_with_GS-rich_pct": f"{100 * len(cls_prot['GS-rich']) / n_prot:.1f}",
        })
        print(f"[done] {key}: {len(rows):,} regions in {len(hits):,} proteins")
    P.write_tsv(os.path.join(a.out, "lcr_summary.tsv"), summary, SUMMARY_FIELDS)
    return 0


if __name__ == "__main__":
    sys.exit(main())
