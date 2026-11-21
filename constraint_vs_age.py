#!/usr/bin/env python3
"""Step 7a — constraint vs tract age.

Joins the dated P. falciparum poly-N tracts (step 6) with gene-level
essentiality and expression:
  - P. falciparum piggyBac (Zhang et al. 2018): Mutagenesis Index Score
    (MIS; low = essential; Zhang's essential call MIS <= 0.26 used here as
    'essential') and Mutant Fitness Score (MFS).
  - P. berghei PlasmoGEM (Bushell et al. 2017): phenotype of the P. berghei
    reciprocal-best-hit ortholog (Essential / Slow / Dispensable).
  - Expression: maximum RNA-seq percentile across asexual stages (as in the
    Paper 3 annotation table).

Questions
  Q1 gene level: are genes carrying old tracts (present in all Laverania) more
     often essential / more highly expressed than genes whose tracts are all
     young (Pf only or Pf+Pprf)? Background: genes with no N-rich tract.
  Q2 tract level: within essential vs dispensable genes, is the age
     distribution of tracts different? (retention under constraint would make
     tracts in essential genes older; purging would make them rarer/younger)
  Q3 deep retention: tracts whose site is also N-rich outside Laverania vs
     Laverania-only tracts, in essential vs dispensable genes.
Protein length is a confounder (long proteins carry more tracts and differ in
MIS), so Q1 is repeated within length quartiles and with a logistic model.
"""
import argparse
import csv
import json
import os
from collections import Counter, defaultdict

import numpy as np
from scipy import stats

YOUNG = {"Pf only", "Pf+Pprf"}
OLD = {"all Laverania (+clade A)"}
MID = {"+Pr", "+Pbilc", "clade B (+Pblac)"}


def rd(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def fisher(a, b, c, d):
    o, p = stats.fisher_exact([[a, b], [c, d]])
    return round(float(o), 3) if np.isfinite(o) else None, float(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dating", required=True)
    ap.add_argument("--pf-annot", required=True)
    ap.add_argument("--pb-annot", required=True)
    ap.add_argument("--expression", required=True)
    ap.add_argument("--rbh-pb", required=True, help="rbh_pfalciparum__pberghei.m8")
    ap.add_argument("--pf-proteins", required=True)
    ap.add_argument("--mis-essential", type=float, default=0.26)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import polylib as P

    recs = [(P.parse_veupath_header(h), s.rstrip("*")) for h, s in P.read_fasta(a.pf_proteins)]
    prim = P.primary_isoforms([(m, s) for m, s in recs if not P.is_pseudo(m)])
    plen = {g: len(s) for g, (_, s) in prim.items()}

    mis, mfs = {}, {}
    for r in rd(a.pf_annot):
        g = r["Gene ID"]
        mis[g] = num(r.get("Mutagenesis Index Score"))
        mfs[g] = num(r.get("Mutant Fitness Score"))
    pbph = {r["Gene ID"]: r.get("Phenotype") for r in rd(a.pb_annot)
            if r.get("Phenotype") not in (None, "", "N/A")}
    expr = {r["gene"]: num(r["value"]) for r in rd(a.expression) if r["dataset"] == "pf"}
    pf2pb = {}
    with open(a.rbh_pb) as fh:
        for line in fh:
            f = line.split("\t")
            pf2pb[f[0]] = f[1]

    tracts = [r for r in rd(a.dating) if r["origin"] != "uninformative"]
    by_gene = defaultdict(list)
    for t in tracts:
        by_gene[t["gene"]].append(t)

    def gene_class(ts):
        o = {t["origin"] for t in ts}
        if o & OLD and not o & YOUNG:
            return "old_only"
        if o <= YOUNG:
            return "young_only"
        return "mixed"

    genes = [g for g in plen if mis.get(g) is not None]
    gclass = {g: (gene_class(by_gene[g]) if g in by_gene else "no_tract") for g in genes}

    def summarise(gs):
        m = [mis[g] for g in gs]
        e = [expr[g] for g in gs if expr.get(g) is not None]
        pb = [pbph.get(pf2pb.get(g)) for g in gs if pbph.get(pf2pb.get(g))]
        return {"n": len(gs), "pct_essential_MIS": round(100 * np.mean([x <= a.mis_essential for x in m]), 1),
                "median_MIS": round(float(np.median(m)), 3),
                "median_MFS": round(float(np.median([mfs[g] for g in gs if mfs.get(g) is not None])), 3),
                "median_expr_pct": round(float(np.median(e)), 1) if e else None,
                "pb_n": len(pb), "pb_pct_essential": round(100 * np.mean([x == "Essential" for x in pb]), 1) if pb else None,
                "median_protein_len": int(np.median([plen[g] for g in gs]))}

    res = {"Q1_gene_classes": {c: summarise([g for g in genes if gclass[g] == c])
                               for c in ("no_tract", "young_only", "mixed", "old_only")}}

    # Q1 tests
    old = [g for g in genes if gclass[g] == "old_only"]
    yng = [g for g in genes if gclass[g] == "young_only"]
    none = [g for g in genes if gclass[g] == "no_tract"]
    ess = lambda gs: sum(mis[g] <= a.mis_essential for g in gs)
    res["Q1_tests"] = {
        "old_vs_young_essential_fisher": fisher(ess(old), len(old) - ess(old), ess(yng), len(yng) - ess(yng)),
        "old_vs_young_MIS_mannwhitney_p": float(stats.mannwhitneyu([mis[g] for g in old], [mis[g] for g in yng]).pvalue),
        "tract_vs_none_essential_fisher": fisher(ess(old + yng), len(old + yng) - ess(old + yng), ess(none), len(none) - ess(none)),
        "old_vs_young_expr_mannwhitney_p": float(stats.mannwhitneyu(
            [expr[g] for g in old if expr.get(g) is not None],
            [expr[g] for g in yng if expr.get(g) is not None]).pvalue),
    }

    # length-stratified essential fraction (quartiles over all genes)
    qs = np.quantile([plen[g] for g in genes], [0.25, 0.5, 0.75])
    def lq(g):
        return int(np.searchsorted(qs, plen[g]))
    strat = {}
    for q in range(4):
        row = {}
        for c in ("no_tract", "young_only", "mixed", "old_only"):
            gs = [g for g in genes if gclass[g] == c and lq(g) == q]
            row[c] = {"n": len(gs), "pct_essential": round(100 * np.mean([mis[g] <= a.mis_essential for g in gs]), 1) if gs else None}
        strat[f"length_Q{q + 1}"] = row
    res["Q1_length_stratified"] = strat

    # logistic: essential ~ log length + class dummies (old_only, young_only, mixed; ref no_tract)
    try:
        import statsmodels.api as sm
        X, y = [], []
        for g in genes:
            c = gclass[g]
            X.append([np.log(plen[g]), c == "young_only", c == "mixed", c == "old_only"])
            y.append(mis[g] <= a.mis_essential)
        X = sm.add_constant(np.array(X, dtype=float))
        fit = sm.Logit(np.array(y, dtype=float), X).fit(disp=0)
        names = ["const", "log_length", "young_only", "mixed", "old_only"]
        res["Q1_logistic_essential"] = {n: {"coef": round(float(c), 3), "OR": round(float(np.exp(c)), 3),
                                            "p": float(p)} for n, c, p in zip(names, fit.params, fit.pvalues)}
    except Exception as exc:  # noqa: BLE001
        res["Q1_logistic_essential"] = f"statsmodels unavailable: {exc}"

    # Q2 tract level: age distribution in essential vs dispensable genes
    def age_dist(ts):
        c = Counter("young" if t["origin"] in YOUNG else "old" if t["origin"] in OLD else "mid" for t in ts)
        n = sum(c.values())
        return {k: round(100 * c[k] / n, 1) for k in ("young", "mid", "old")} | {"n": n}
    t_ess = [t for t in tracts if mis.get(t["gene"]) is not None and mis[t["gene"]] <= a.mis_essential]
    t_dis = [t for t in tracts if mis.get(t["gene"]) is not None and mis[t["gene"]] > 0.8]
    res["Q2_tract_age_by_essentiality"] = {"essential_MIS<=0.26": age_dist(t_ess),
                                           "dispensable_MIS>0.8": age_dist(t_dis)}
    yo = lambda ts: sum(t["origin"] in YOUNG for t in ts)
    res["Q2_young_fraction_fisher"] = fisher(yo(t_ess), len(t_ess) - yo(t_ess), yo(t_dis), len(t_dis) - yo(t_dis))

    # Q3 deep retention outside Laverania
    deep = lambda t: bool(t.get("outside_laverania"))
    for lbl, ts in (("essential", t_ess), ("dispensable", t_dis)):
        o = [t for t in ts if t["origin"] in OLD]
        res.setdefault("Q3_old_tracts_shared_outside_laverania", {})[lbl] = {
            "n_old": len(o), "pct_shared_outside": round(100 * np.mean([deep(t) for t in o]), 1) if o else None}
    oe = [t for t in t_ess if t["origin"] in OLD]
    od = [t for t in t_dis if t["origin"] in OLD]
    res["Q3_fisher"] = fisher(sum(map(deep, oe)), len(oe) - sum(map(deep, oe)),
                              sum(map(deep, od)), len(od) - sum(map(deep, od)))

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "constraint_vs_age.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
