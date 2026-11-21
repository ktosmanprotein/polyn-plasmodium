#!/usr/bin/env python3
"""Step 1b — catalogue poly-N runs and low-complexity regions, and compare
their abundance with two null models.

Per species this writes
  <key>_proteins.tsv   one row per protein (primary isoform, non-pseudo)
  <key>_tracts.tsv     one row per tract (poly-N run or low-complexity region)
and across species
  species_summary.tsv  observed vs expected run counts under both nulls
  SUMMARY.md           the same, readable

Null models
  composition  Each protein keeps its own residue frequencies; residues are
               placed independently. Asks: are N residues clustered into runs
               beyond what this protein's N content alone predicts?
  mutation     Every residue is drawn from codons built from the species'
               nucleotide composition at fourfold-degenerate sites (the usual
               proxy for mutational equilibrium). Asks: is there more poly-N
               than the genome's AT bias alone would produce?
Lysine is run through both nulls as a control: its codons (AAA/AAG) are as
A-rich as asparagine's (AAT/AAC), so mutation bias predicts similar amounts.
"""
import argparse
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polylib as P  # noqa: E402

PROTEIN_FIELDS = [
    "gene", "uniprot", "length", "n_frac", "k_frac", "longest_n_run",
    "n_runs_short", "n_runs_long", "k_runs_short", "k_runs_long",
    "lcr_n", "lcr_residues", "lcr_frac", "lcr_classes", "has_polyN",
]
TRACT_FIELDS = [
    "gene", "uniprot", "tract_type", "tract_class", "start", "end", "length",
    "n_frac", "composition", "rel_position",
]
SUMMARY_FIELDS = [
    "species", "label", "proteins", "mapped_uniprot", "aa_residues",
    "fourfold_sites", "fourfold_AT",
    "N_freq_obs", "N_freq_mut", "K_freq_obs", "K_freq_mut",
    "N_runs_short_obs", "N_runs_short_exp_comp", "N_runs_short_exp_mut",
    "N_runs_long_obs", "N_runs_long_exp_comp", "N_runs_long_exp_mut",
    "K_runs_short_obs", "K_runs_short_exp_comp", "K_runs_short_exp_mut",
    "K_runs_long_obs", "K_runs_long_exp_comp", "K_runs_long_exp_mut",
    "proteins_with_polyN", "proteins_with_polyN_pct",
    "proteins_with_LCR_pct", "lcr_residue_pct",
    "lcr_N-rich", "lcr_Q-rich", "lcr_GS-rich", "lcr_charged-KED", "lcr_mixed",
]


def uniprot_map(apico_dir, key):
    """gene id -> UniProt accession, from the apicofold proteome table."""
    path = os.path.join(apico_dir, f"{key}.tsv") if apico_dir else ""
    out = {}
    if not path or not os.path.exists(path):
        return out
    for r in P.read_tsv(path):
        for tok in (r.get("veupathdb") or "").split(";"):
            tok = tok.strip()
            if not tok:
                continue
            gene = tok.split(":", 1)[-1]
            out.setdefault(gene, r["accession"])
    return out


def load_proteins(path):
    recs = []
    dropped = Counter()
    for h, s in P.read_fasta(path):
        meta = P.parse_veupath_header(h)
        if P.is_pseudo(meta):
            dropped["pseudo"] += 1
            continue
        s = s.upper().rstrip("*")
        if "*" in s:
            dropped["internal_stop"] += 1
            continue
        if not s:
            dropped["empty"] += 1
            continue
        recs.append((meta, s))
    return P.primary_isoforms(recs), dropped


def load_cds(path, keep_transcripts):
    for h, s in P.read_fasta(path):
        meta = P.parse_veupath_header(h)
        if P.is_pseudo(meta):
            continue
        if keep_transcripts and meta["transcript"] not in keep_transcripts \
                and meta["id"] not in keep_transcripts:
            continue
        yield s


def analyse_species(sp, args):
    key = sp["key"]
    prot_path = os.path.join(args.seqs, f"{key}_proteins.fasta")
    cds_path = os.path.join(args.seqs, f"{key}_cds.fasta")
    if not os.path.exists(prot_path):
        print(f"[skip] {key}: {prot_path} missing (run step 01)")
        return None

    proteins, dropped = load_proteins(prot_path)
    to_uniprot = uniprot_map(args.apico_proteomes, key)

    # Mutation-bias null from fourfold sites of the same genes' CDS.
    if os.path.exists(cds_path):
        keep = {m["transcript"] for m, _ in proteins.values()}
        four, n_sites = P.fourfold_composition(load_cds(cds_path, keep))
    else:
        print(f"[warn] {key}: no CDS file; mutation null unavailable")
        four, n_sites = None, 0
    p_mut = {aa: (P.residue_prob_under_composition(four, P.CODONS[aa])
                  if four else float("nan")) for aa in ("N", "K")}

    rules = P.DEFAULT_CLASS_RULES
    ks, kl = args.run_short, args.run_long
    prot_rows, tract_rows = [], []
    tot = Counter()
    exp = defaultdict(float)
    lcr_class_counts = Counter()

    for gene, (meta, seq) in sorted(proteins.items()):
        L = len(seq)
        comp = Counter(seq)
        pN, pK = comp["N"] / L, comp["K"] / L
        uni = to_uniprot.get(gene, "")

        runs = {}
        for aa in ("N", "K"):
            runs[aa] = P.homopolymer_runs(seq, aa, ks)
            long_n = sum(1 for s, e in runs[aa] if e - s >= kl)
            tot[f"{aa}_short"] += len(runs[aa])
            tot[f"{aa}_long"] += long_n
            pa = pN if aa == "N" else pK
            exp[f"{aa}_short_comp"] += P.expected_maximal_runs(L, pa, ks)
            exp[f"{aa}_long_comp"] += P.expected_maximal_runs(L, pa, kl)
            exp[f"{aa}_short_mut"] += P.expected_maximal_runs(L, p_mut[aa], ks)
            exp[f"{aa}_long_mut"] += P.expected_maximal_runs(L, p_mut[aa], kl)

        lcrs = P.lcr_regions(seq, args.lcr_window, args.lcr_trigger,
                             args.lcr_extend, args.lcr_min_len)
        lcr_res = sum(e - s for s, e in lcrs)
        classes = []
        for s, e in lcrs:
            seg = seq[s:e]
            cls = P.classify_segment(seg, rules)
            classes.append(cls)
            lcr_class_counts[cls] += 1
            tract_rows.append({
                "gene": gene, "uniprot": uni, "tract_type": "LCR",
                "tract_class": cls, "start": s + 1, "end": e, "length": e - s,
                "n_frac": f"{seg.count('N') / len(seg):.3f}",
                "composition": P.top_composition(seg),
                "rel_position": f"{((s + e) / 2) / L:.3f}",
            })
        for s, e in runs["N"]:
            tract_rows.append({
                "gene": gene, "uniprot": uni, "tract_type": "polyN",
                "tract_class": "N-run", "start": s + 1, "end": e,
                "length": e - s, "n_frac": "1.000", "composition": "N:1.00",
                "rel_position": f"{((s + e) / 2) / L:.3f}",
            })

        tot["proteins"] += 1
        tot["residues"] += L
        tot["N"] += comp["N"]
        tot["K"] += comp["K"]
        tot["mapped"] += bool(uni)
        tot["with_polyN"] += bool(runs["N"])
        tot["with_lcr"] += bool(lcrs)
        tot["lcr_residues"] += lcr_res
        prot_rows.append({
            "gene": gene, "uniprot": uni, "length": L,
            "n_frac": f"{pN:.4f}", "k_frac": f"{pK:.4f}",
            "longest_n_run": max_run(seq, "N"),
            "n_runs_short": len(runs["N"]),
            "n_runs_long": sum(1 for s, e in runs["N"] if e - s >= kl),
            "k_runs_short": len(runs["K"]),
            "k_runs_long": sum(1 for s, e in runs["K"] if e - s >= kl),
            "lcr_n": len(lcrs), "lcr_residues": lcr_res,
            "lcr_frac": f"{lcr_res / L:.4f}",
            "lcr_classes": ",".join(classes),
            "has_polyN": int(bool(runs["N"])),
        })

    P.write_tsv(os.path.join(args.out, f"{key}_proteins.tsv"), prot_rows, PROTEIN_FIELDS)
    P.write_tsv(os.path.join(args.out, f"{key}_tracts.tsv"), tract_rows, TRACT_FIELDS)

    n = tot["proteins"]
    res = tot["residues"]
    row = {
        "species": key, "label": sp.get("label", key), "proteins": n,
        "mapped_uniprot": tot["mapped"], "aa_residues": res,
        "fourfold_sites": n_sites,
        "fourfold_AT": f"{(four['A'] + four['T']):.4f}" if four else "",
        "N_freq_obs": f"{tot['N'] / res:.4f}", "N_freq_mut": f"{p_mut['N']:.4f}",
        "K_freq_obs": f"{tot['K'] / res:.4f}", "K_freq_mut": f"{p_mut['K']:.4f}",
        "proteins_with_polyN": tot["with_polyN"],
        "proteins_with_polyN_pct": f"{100 * tot['with_polyN'] / n:.1f}",
        "proteins_with_LCR_pct": f"{100 * tot['with_lcr'] / n:.1f}",
        "lcr_residue_pct": f"{100 * tot['lcr_residues'] / res:.1f}",
    }
    for aa in ("N", "K"):
        for size in ("short", "long"):
            row[f"{aa}_runs_{size}_obs"] = tot[f"{aa}_{size}"]
            row[f"{aa}_runs_{size}_exp_comp"] = f"{exp[f'{aa}_{size}_comp']:.4g}"
            row[f"{aa}_runs_{size}_exp_mut"] = f"{exp[f'{aa}_{size}_mut']:.4g}"
    for cls in ("N-rich", "Q-rich", "GS-rich", "charged-KED", "mixed"):
        row[f"lcr_{cls}"] = lcr_class_counts[cls]

    drop_txt = ", ".join(f"{v} {k}" for k, v in dropped.items()) or "none"
    print(f"[done] {key}: {n:,} proteins ({tot['mapped']:,} mapped to UniProt; "
          f"dropped: {drop_txt}); {tot['with_polyN']:,} with an N run >= {ks}")
    return row


def max_run(seq, aa):
    best = cur = 0
    for c in seq:
        cur = cur + 1 if c == aa else 0
        if cur > best:
            best = cur
    return best


def ratio(obs, exp):
    try:
        exp = float(exp)
        if exp <= 0:
            return "inf"
        r = float(obs) / exp
        return f"{r:,.1f}x" if r < 1e4 else f"{r:.1e}x"
    except ValueError:
        return "n/a"


def write_markdown(path, rows, args):
    ks, kl = args.run_short, args.run_long
    lines = [
        "# Step 1 — poly-N runs and low-complexity regions",
        "",
        f"Runs: maximal homopolymers of length >= {ks} (short) and >= {kl} (long).",
        f"LCRs: SEG-like, window {args.lcr_window}, trigger {args.lcr_trigger} bits, "
        f"extend {args.lcr_extend} bits, min length {args.lcr_min_len}.",
        "",
        "## Observed vs expected runs",
        "",
        "obs/exp under the **composition** null (each protein's own residue "
        "frequencies) and the **mutation** null (fourfold-site nucleotide "
        "composition). Lysine is the control.",
        "",
        f"| Species | 4-fold AT | N freq obs / mut | N runs ≥{ks} | ÷comp | ÷mut | "
        f"N runs ≥{kl} | ÷comp | ÷mut | K runs ≥{ks} | ÷comp | ÷mut |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['label']} | {r['fourfold_AT']} | {r['N_freq_obs']} / {r['N_freq_mut']} | "
            f"{r['N_runs_short_obs']:,} | {ratio(r['N_runs_short_obs'], r['N_runs_short_exp_comp'])} | "
            f"{ratio(r['N_runs_short_obs'], r['N_runs_short_exp_mut'])} | "
            f"{r['N_runs_long_obs']:,} | {ratio(r['N_runs_long_obs'], r['N_runs_long_exp_comp'])} | "
            f"{ratio(r['N_runs_long_obs'], r['N_runs_long_exp_mut'])} | "
            f"{r['K_runs_short_obs']:,} | {ratio(r['K_runs_short_obs'], r['K_runs_short_exp_comp'])} | "
            f"{ratio(r['K_runs_short_obs'], r['K_runs_short_exp_mut'])} |")
    lines += [
        "",
        "## Tract burden",
        "",
        f"| Species | proteins | with N run ≥{ks} | with LCR | LCR residues | "
        "N-rich | Q-rich | GS-rich | charged | mixed |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['label']} | {r['proteins']:,} | {r['proteins_with_polyN']:,} "
            f"({r['proteins_with_polyN_pct']}%) | {r['proteins_with_LCR_pct']}% | "
            f"{r['lcr_residue_pct']}% | {r['lcr_N-rich']:,} | {r['lcr_Q-rich']:,} | "
            f"{r['lcr_GS-rich']:,} | {r['lcr_charged-KED']:,} | {r['lcr_mixed']:,} |")
    lines += [
        "",
        "## Reading the nulls",
        "",
        "- ÷comp well above 1: N residues are clustered into runs beyond what "
        "each protein's N content predicts (runs are not scattered Asn).",
        "- ÷mut well above 1 for N but near 1 for K: the excess is specific to "
        "asparagine, not a by-product of AT-rich codons.",
        "- The mutation null treats residues as independent; it is a floor for "
        "what AT bias alone produces, not a full model of protein evolution.",
        "",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--seqs", required=True)
    ap.add_argument("--apico-proteomes", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--run-short", type=int, default=5)
    ap.add_argument("--run-long", type=int, default=10)
    ap.add_argument("--lcr-window", type=int, default=12)
    ap.add_argument("--lcr-trigger", type=float, default=2.2)
    ap.add_argument("--lcr-extend", type=float, default=2.5)
    ap.add_argument("--lcr-min-len", type=int, default=15)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    rows = []
    for sp in P.read_tsv(args.species):
        r = analyse_species(sp, args)
        if r:
            rows.append(r)
    if not rows:
        print("Nothing analysed.")
        return 1
    P.write_tsv(os.path.join(args.out, "species_summary.tsv"), rows, SUMMARY_FIELDS)
    md = os.path.join(args.out, "SUMMARY.md")
    write_markdown(md, rows, args)
    print(f"\nSummary: {md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
