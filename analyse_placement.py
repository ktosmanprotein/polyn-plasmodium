#!/usr/bin/env python3
"""Step 5 — tests on the placement tables from step 4.

A. Where disordered low-complexity regions sit (location by class).
B. Composition turnover vs conservation at homologous disordered positions.
   For every disordered region in A whose ortholog B has its own disordered
   segment at the homologous position (internal_variable / terminal_variable),
   ask whether B's segment carries a region of the SAME class. Null: B
   segments permuted among A regions within the same ordered species pair and
   the same B-segment length bin (longer segments are more likely to contain
   any region), 2,000 permutations. A continuous version uses the asparagine
   fraction of B's segment for A regions that are N-rich.
C. Lineage pattern of disordered LCR presence across the four Plasmodium
   orthologs, observed vs expected from each species' own occupancy rate.
   With one rodent species, fever exposure and the P. berghei lineage are
   confounded; this is descriptive.
D. Robustness: B repeated on orthogroups whose models are all < 1,500 aa,
   because AlphaFold DB lacks many long proteins, unevenly by species.

Insertion-site geometry (internal_insertion) needs the structures and is in
insertion_site_null.py.
"""
import argparse
import csv
import json
import os
from collections import Counter, defaultdict

import numpy as np

PLAS = ["pfalciparum", "pvivax", "pknowlesi", "pberghei"]
SHORT = {"pfalciparum": "Pf", "pvivax": "Pv", "pknowlesi": "Pk", "pberghei": "Pb", "tgondii": "Tg"}
CLASSES = ["N-rich", "GS-rich", "charged", "Q-rich", "other"]
HOMOL = ("internal_variable", "terminal_variable")


def rd(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def length_bin(n, edges=(20, 40, 80, 160, 320)):
    return sum(n > e for e in edges)


def permutation_test(rows, score_fn, match_fn, n_perm, rng):
    """Null for 'does B at the homologous site resemble A?'.

    rows: all homologous-site rows (every A class). Within each stratum (same
    ordered species pair, same B-segment length bin) the B sides are shuffled
    across A regions; the statistic is the mean of match_fn(A, B) over the A
    rows selected by score_fn. Returns (observed, null array)."""
    strata = defaultdict(list)
    for r in rows:
        strata[(r["a_species"], r["b_species"], length_bin(int(r["b_gap"] or 0)))].append(r)
    groups = list(strata.values())
    sel = [[i for i, r in enumerate(g) if score_fn(r)] for g in groups]
    n_sel = sum(len(x) for x in sel)
    if n_sel == 0:
        return float("nan"), np.array([np.nan])
    obs = sum(match_fn(g[i], g[i]) for g, ix in zip(groups, sel) for i in ix) / n_sel
    null = np.empty(n_perm)
    for k in range(n_perm):
        tot = 0.0
        for g, ix in zip(groups, sel):
            if not ix:
                continue
            perm = rng.permutation(len(g))
            tot += sum(match_fn(g[i], g[perm[i]]) for i in ix)
        null[k] = tot / n_sel
    return obs, null


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--placement", required=True, help="results/placement_v3")
    ap.add_argument("--out", required=True)
    ap.add_argument("--perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--disorder-plddt", type=float, default=70.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    rp = rd(os.path.join(a.placement, "region_partner.tsv"))
    mem = rd(os.path.join(a.placement, "members.tsv"))
    out = {}

    pp = [r for r in rp if r["a_species"] in PLAS and r["b_species"] in PLAS]
    dis = [r for r in pp if float(r["region_plddt"]) < a.disorder_plddt]

    # ---- A. location by class --------------------------------------------
    locA = {}
    for cls in CLASSES:
        c = Counter(r["location"] for r in dis if r["lcr_class"] == cls)
        locA[cls] = dict(c)
    out["A_location_by_class"] = locA
    out["A_n_pairs_plasmodium"] = len(pp)
    out["A_n_pairs_disordered"] = len(dis)
    out["A_poor_alignment_frac"] = sum(r["location"] == "poor_alignment" for r in dis) / len(dis)

    # ---- B. turnover vs conservation -------------------------------------
    hom = [r for r in dis if r["location"] in HOMOL and r["b_seg_class"] != ""]

    def same_class(ar, br):
        return 1.0 if br["b_seg_class"] == ar["lcr_class"] else 0.0

    def any_lcr(ar, br):
        return 1.0 if br["b_seg_class"] not in ("none", "") else 0.0

    resB = {}
    for cls in CLASSES:
        n_cls = sum(r["lcr_class"] == cls for r in hom)
        if n_cls < 30:
            continue
        o, nl = permutation_test(hom, lambda r, c=cls: r["lcr_class"] == c, same_class, a.perm, rng)
        resB[cls] = {"n": n_cls, "observed_same_class": round(float(o), 4),
                     "null_mean": round(float(nl.mean()), 4),
                     "null_95": [round(float(x), 4) for x in np.percentile(nl, [2.5, 97.5])],
                     "ratio": round(float(o / nl.mean()), 2) if nl.mean() > 0 else None,
                     "p_one_sided": float((np.sum(nl >= o) + 1) / (len(nl) + 1))}
    out["B_same_class_at_homologous_site"] = resB
    o, nl = permutation_test(hom, lambda r: True, any_lcr, a.perm, rng)
    out["B_any_lcr_at_homologous_site"] = {"n": len(hom), "observed": round(float(o), 4),
                                           "null_mean": round(float(nl.mean()), 4)}

    # continuous: N fraction of B segment when A region is N-rich
    hom_n = [r for r in hom if r["b_seg_n_frac"]]
    nr = [r for r in hom_n if r["lcr_class"] == "N-rich"]
    o, nl = permutation_test(hom_n, lambda r: r["lcr_class"] == "N-rich",
                             lambda ar, br: float(br["b_seg_n_frac"]), a.perm, rng)
    allb = [float(r["b_seg_n_frac"]) for r in hom if r["b_seg_n_frac"]]
    out["B_Nfrac_partner_segment_when_A_Nrich"] = {
        "n": len(nr), "observed_mean": round(float(o), 4), "null_mean": round(float(nl.mean()), 4),
        "p_one_sided": float((np.sum(nl >= o) + 1) / (len(nl) + 1)),
        "all_homologous_segments_mean": round(float(np.mean(allb)), 4)}

    # species-pair matrix for N-rich: fraction of A N-rich regions whose
    # homologous B segment is also N-rich
    mat = {}
    for sa in PLAS:
        for sb in PLAS:
            if sa == sb:
                continue
            rows = [r for r in hom if r["lcr_class"] == "N-rich" and r["a_species"] == sa
                    and r["b_species"] == sb]
            if rows:
                c = Counter(r["b_seg_class"] for r in rows)
                mat[f"{SHORT[sa]}->{SHORT[sb]}"] = {"n": len(rows),
                                                   **{k: round(v / len(rows), 3) for k, v in c.items()}}
    out["B_Nrich_partner_class_matrix"] = mat

    # ---- C. lineage pattern of disordered LCR presence --------------------
    has = defaultdict(lambda: defaultdict(set))   # og -> species -> classes
    placed = defaultdict(set)
    for r in mem:
        if r["species"] in PLAS:
            placed[r["og_id"]].add(r["species"])
    regs = rd(os.path.join(a.placement, "regions.tsv"))
    for r in regs:
        if r["species"] in PLAS and float(r["mean_plddt"]) < a.disorder_plddt:
            has[r["og_id"]][r["species"]].add(r["lcr_class"])
    full = [og for og, s in placed.items() if len(s) == 4]
    resC = {}
    for label, test in (("any", lambda cl: bool(cl)), ("N-rich", lambda cl: "N-rich" in cl),
                        ("GS-rich", lambda cl: "GS-rich" in cl)):
        pats = Counter()
        for og in full:
            pats["".join("1" if test(has[og][s]) else "0" for s in PLAS)] += 1
        n = len(full)
        p = {s: sum(int(k[i]) * v for k, v in pats.items()) / n for i, s in enumerate(PLAS)}
        rows = []
        for k in sorted(pats, key=lambda x: -pats[x]):
            e = n
            for i, s in enumerate(PLAS):
                e *= p[s] if k[i] == "1" else 1 - p[s]
            rows.append({"pattern": k, "observed": pats[k], "expected_independent": round(e, 1),
                         "o_over_e": round(pats[k] / e, 2) if e else None})
        resC[label] = {"n_orthogroups": n, "occupancy": {SHORT[s]: round(v, 3) for s, v in p.items()},
                       "patterns": rows}
    out["C_lineage_patterns"] = resC
    out["C_pattern_order"] = [SHORT[s] for s in PLAS]

    # ---- D. robustness: drop orthogroups with any model >= 1500 aa --------
    long_og = {r["og_id"] for r in mem if int(r["model_len"]) >= 1500}
    homD = [r for r in hom if r["og_id"] not in long_og]
    resD = {}
    for cls in CLASSES:
        n_cls = sum(r["lcr_class"] == cls for r in homD)
        if n_cls < 30:
            continue
        o, nl = permutation_test(homD, lambda r, c=cls: r["lcr_class"] == c, same_class,
                                 max(500, a.perm // 4), rng)
        resD[cls] = {"n": n_cls, "observed_same_class": round(float(o), 4),
                     "null_mean": round(float(nl.mean()), 4),
                     "ratio": round(float(o / nl.mean()), 2) if nl.mean() > 0 else None}
    locD = {cls: dict(Counter(r["location"] for r in dis if r["lcr_class"] == cls
                              and r["og_id"] not in long_og)) for cls in CLASSES}
    out["D_short_proteins_only"] = {"orthogroups_removed": len(long_og),
                                    "same_class": resD, "location_by_class": locD}

    with open(os.path.join(a.out, "step5_results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1)[:6000])


if __name__ == "__main__":
    main()
