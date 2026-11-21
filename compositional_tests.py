#!/usr/bin/env python3
"""Step 5d — genome-level vs protein-level tests on poly-N.

Test 1  Asparagine vs lysine against mutation bias.
        AAT/AAC (N) and AAA/AAG (K) share their first two bases, so under
        mutation bias alone the expected N:K ratio is (T+C)/(A+G) at the
        third codon position, estimated at fourfold-degenerate sites of the
        same species. Observed N:K is measured in: the whole proteome; the
        disordered segments that sit at homologous positions (step 4); those
        segments with the low-complexity tracts removed (the 'background'
        disorder in which tracts arise); the tracts themselves.
        A site-matched version compares, at the same homologous segments,
        P. falciparum with each other species.
Test 2  Seeding. For every P. falciparum N-rich tract, what does the
        ortholog's segment at the same position contain (N and K content,
        short N runs, AAT/AAC and AAA/AAG codon runs, its own region class)?
        Control: ortholog segments opposite P. falciparum tracts of other
        classes, same partner species.
Test 4  Asymmetry against a background null. Fraction of A's N-rich tracts
        whose homologous B segment is N-rich, for every ordered species pair,
        with B segments permuted among homologous sites of the same species
        pair and length bin (B's base rate is kept). O/E > 1 means real
        positional association beyond B's overall poly-N frequency.

Coordinates: only proteins whose AlphaFold model sequence is identical to
the annotation (step 4 members.tsv, seq_identity == 1) are used, so model and
annotation positions coincide and CDS codons can be read directly.
"""
import argparse
import csv
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polylib as P  # noqa: E402
from analyse_placement import permutation_test  # noqa: E402

PLAS = ["pfalciparum", "pvivax", "pknowlesi", "pberghei"]
SHORT = {"pfalciparum": "Pf", "pvivax": "Pv", "pknowlesi": "Pk", "pberghei": "Pb"}
HOMOL = ("internal_variable", "terminal_variable")


def rd(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_species(seqs_dir, sp):
    recs = []
    for h, s in P.read_fasta(os.path.join(seqs_dir, f"{sp}_proteins.fasta")):
        m = P.parse_veupath_header(h)
        if not P.is_pseudo(m):
            recs.append((m, s.upper().rstrip("*")))
    prim = P.primary_isoforms(recs)
    tx = {m["transcript"]: g for g, (m, _) in prim.items()}
    cds = {}
    for h, s in P.read_fasta(os.path.join(seqs_dir, f"{sp}_cds.fasta")):
        m = P.parse_veupath_header(h)
        g = tx.get(m["transcript"]) or tx.get(m["id"])
        if g:
            cds[g] = s.upper()
    prot = {g: s for g, (_, s) in prim.items()}
    # keep CDS only where it translates to the protein length (+stop)
    good = {g: c for g, c in cds.items() if len(c) in (3 * len(prot[g]), 3 * len(prot[g]) + 3)}
    return prot, good


def third_position_comp(cds_dict):
    f, n = P.fourfold_composition(cds_dict.values())
    return f, n


def nk(seq):
    return seq.count("N"), seq.count("K")


def longest_codon_run(cds, codons):
    best = cur = 0
    for i in range(0, len(cds) - 2, 3):
        if cds[i:i + 3] in codons:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def longest_run(seq, aa):
    m = re.findall(f"{aa}+", seq)
    return max((len(x) for x in m), default=0)


def a_segment(r, LA):
    """A's homologous segment, 0-based [s, e)."""
    if r["location"] == "internal_variable":
        return int(r["a_anchor_left"]), int(r["a_anchor_right"]) - 1
    if r["a_anchor_right"] and not r["a_anchor_left"]:
        return 0, int(r["a_anchor_right"]) - 1
    return int(r["a_anchor_left"]), LA


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--placement", required=True)
    ap.add_argument("--seqs", required=True)
    ap.add_argument("--lcr-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    out = {}

    prot, cds, comp = {}, {}, {}
    for sp in PLAS:
        prot[sp], cds[sp] = load_species(a.seqs, sp)
        f, n = third_position_comp(cds[sp])
        comp[sp] = {"fourfold_freq": {b: round(f[b], 4) for b in "ACGT"}, "sites": n,
                    "expected_N_over_K": round((f["T"] + f["C"]) / (f["A"] + f["G"]), 4)}
    out["fourfold_composition"] = comp

    mem = {(r["species"], r["gene"]): r for r in rd(os.path.join(a.placement, "members.tsv"))}

    def exact(sp, g):
        m = mem.get((sp, g))
        return (m is not None and float(m["seq_identity"]) >= 0.9999 and g in prot[sp]
                and int(m["model_len"]) == len(prot[sp][g]))

    rp = [r for r in rd(os.path.join(a.placement, "region_partner.tsv"))
          if r["a_species"] in PLAS and r["b_species"] in PLAS
          and float(r["region_plddt"]) < 70]
    hom = [r for r in rp if r["location"] in HOMOL and r["b_seg_start"]
           and exact(r["a_species"], r["a_gene"]) and exact(r["b_species"], r["b_gene"])]
    out["n_homologous_rows_exact_coords"] = len(hom)

    lcr = defaultdict(list)
    for sp in PLAS:
        for r in rd(os.path.join(a.lcr_dir, f"{sp}_lcr.tsv")):
            lcr[(sp, r["gene"])].append((int(r["start"]) - 1, int(r["end"])))

    # ---- Test 1: compartments ------------------------------------------------
    segs = defaultdict(set)   # species -> {(gene, s, e)}
    for r in hom:
        LA = len(prot[r["a_species"]][r["a_gene"]])
        s, e = a_segment(r, LA)
        if e > s:
            segs[r["a_species"]].add((r["a_gene"], s, e))
        segs[r["b_species"]].add((r["b_gene"], int(r["b_seg_start"]) - 1, int(r["b_seg_end"])))

    t1 = {}
    for sp in PLAS:
        exp = comp[sp]["expected_N_over_K"]
        whole = Counter()
        for s in prot[sp].values():
            whole.update(s)
        seg_n = seg_k = bg_n = bg_k = tr_n = tr_k = 0
        for g, s, e in segs[sp]:
            seq = prot[sp][g]
            mask = [False] * len(seq)
            for x, y in lcr.get((sp, g), []):
                for i in range(max(x, s), min(y, e)):
                    mask[i] = True
            for i in range(s, e):
                c = seq[i]
                if c == "N":
                    seg_n += 1
                    if mask[i]:
                        tr_n += 1
                    else:
                        bg_n += 1
                elif c == "K":
                    seg_k += 1
                    if mask[i]:
                        tr_k += 1
                    else:
                        bg_k += 1

        def row(n, k):
            r = n / k if k else float("inf")
            return {"N": n, "K": k, "N_over_K": round(r, 3),
                    "log2_obs_over_expected": round(math.log2(r / exp), 3) if k and n else None}
        t1[SHORT[sp]] = {"expected_N_over_K": exp,
                         "whole_proteome": row(whole["N"], whole["K"]),
                         "homologous_disordered_segments": row(seg_n, seg_k),
                         "segments_without_tracts": row(bg_n, bg_k),
                         "tracts_in_segments": row(tr_n, tr_k),
                         "n_segments": len(segs[sp])}
    out["test1_N_vs_K"] = t1

    # site-matched: same homologous segments, Pf (A) vs partner (B), tracts removed on both
    def nk_masked(sp, g, s, e):
        seq = prot[sp][g]
        masked = set()
        for x, y in lcr.get((sp, g), []):
            masked.update(range(max(x, s), min(y, e)))
        sub = "".join(seq[i] for i in range(s, e) if i not in masked)
        return sub.count("N"), sub.count("K"), len(sub)

    sm = {}
    for sb in PLAS[1:]:
        dif = []
        for r in hom:
            if r["a_species"] != "pfalciparum" or r["b_species"] != sb:
                continue
            LA = len(prot["pfalciparum"][r["a_gene"]])
            s, e = a_segment(r, LA)
            an, ak, al = nk_masked("pfalciparum", r["a_gene"], s, e)
            bn, bk, bl = nk_masked(sb, r["b_gene"], int(r["b_seg_start"]) - 1, int(r["b_seg_end"]))
            if al >= 10 and bl >= 10:
                dif.append((math.log2((an + 0.5) / (ak + 0.5)), math.log2((bn + 0.5) / (bk + 0.5))))
        d = np.array(dif)
        sm[f"Pf_vs_{SHORT[sb]}"] = {
            "n_sites": len(d),
            "median_log2NK_Pf": round(float(np.median(d[:, 0])), 3),
            "median_log2NK_partner": round(float(np.median(d[:, 1])), 3),
            "median_paired_difference": round(float(np.median(d[:, 0] - d[:, 1])), 3),
            "wilcoxon_p": float(stats.wilcoxon(d[:, 0], d[:, 1]).pvalue),
            "expected_log2NK_Pf": round(math.log2(comp["pfalciparum"]["expected_N_over_K"]), 3),
            "expected_log2NK_partner": round(math.log2(comp[sb]["expected_N_over_K"]), 3)}
    out["test1_site_matched_tracts_removed"] = sm

    # ---- Test 2: seeding ----------------------------------------------------
    def describe(sp, g, s, e):
        seq = prot[sp][g][s:e]
        c = cds[sp].get(g, "")[3 * s:3 * e] if g in cds[sp] else ""
        L = max(len(seq), 1)
        return {"N_frac": seq.count("N") / L, "K_frac": seq.count("K") / L,
                "N_run": longest_run(seq, "N"), "K_run": longest_run(seq, "K"),
                "AAY_run": longest_codon_run(c, {"AAT", "AAC"}) if c else None,
                "AAR_run": longest_codon_run(c, {"AAA", "AAG"}) if c else None}

    t2 = {}
    for sb in PLAS[1:]:
        groups = {"opposite_Pf_Nrich": [], "opposite_Pf_other_classes": []}
        cls = Counter()
        for r in hom:
            if r["a_species"] != "pfalciparum" or r["b_species"] != sb:
                continue
            d = describe(sb, r["b_gene"], int(r["b_seg_start"]) - 1, int(r["b_seg_end"]))
            if r["lcr_class"] == "N-rich":
                groups["opposite_Pf_Nrich"].append(d)
                cls[r["b_seg_class"]] += 1
            else:
                groups["opposite_Pf_other_classes"].append(d)
        res = {}
        for k, v in groups.items():
            res[k] = {"n": len(v),
                      "median_N_frac": round(float(np.median([x["N_frac"] for x in v])), 3),
                      "median_K_frac": round(float(np.median([x["K_frac"] for x in v])), 3),
                      "pct_N_run_ge3": round(100 * np.mean([x["N_run"] >= 3 for x in v]), 1),
                      "pct_K_run_ge3": round(100 * np.mean([x["K_run"] >= 3 for x in v]), 1),
                      "pct_AAY_codon_run_ge3": round(100 * np.mean([(x["AAY_run"] or 0) >= 3 for x in v]), 1),
                      "pct_AAR_codon_run_ge3": round(100 * np.mean([(x["AAR_run"] or 0) >= 3 for x in v]), 1)}
        x = [d["N_frac"] for d in groups["opposite_Pf_Nrich"]]
        y = [d["N_frac"] for d in groups["opposite_Pf_other_classes"]]
        res["N_frac_mannwhitney_p"] = float(stats.mannwhitneyu(x, y, alternative="greater").pvalue)
        xk = [d["K_frac"] for d in groups["opposite_Pf_Nrich"]]
        yk = [d["K_frac"] for d in groups["opposite_Pf_other_classes"]]
        res["K_frac_mannwhitney_p_two_sided"] = float(stats.mannwhitneyu(xk, yk).pvalue)
        res["partner_region_class_opposite_Pf_Nrich"] = {k: round(v / max(sum(cls.values()), 1), 3)
                                                         for k, v in cls.most_common()}
        t2[SHORT[sb]] = res
    out["test2_seeding"] = t2

    # ---- Test 4: asymmetry with background null -----------------------------
    hom_all = [r for r in rp if r["location"] in HOMOL and r["b_seg_class"] != ""]

    def b_is_nrich(ar, br):
        return 1.0 if br["b_seg_class"] == "N-rich" else 0.0
    t4 = {}
    for sa in PLAS:
        for sb in PLAS:
            if sa == sb:
                continue
            rows = [r for r in hom_all if r["a_species"] == sa and r["b_species"] == sb]
            n = sum(r["lcr_class"] == "N-rich" for r in rows)
            if n < 10:
                continue
            o, nl = permutation_test(rows, lambda r: r["lcr_class"] == "N-rich", b_is_nrich,
                                     a.perm, rng)
            t4[f"{SHORT[sa]}->{SHORT[sb]}"] = {
                "n": n, "observed": round(float(o), 3), "expected": round(float(nl.mean()), 3),
                "o_over_e": round(float(o / nl.mean()), 2) if nl.mean() > 0 else None,
                "null_95": [round(float(q), 3) for q in np.percentile(nl, [2.5, 97.5])],
                "p_one_sided": float((np.sum(nl >= o) + 1) / (len(nl) + 1))}
    out["test4_asymmetry_vs_background"] = t4

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "compositional_tests.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
