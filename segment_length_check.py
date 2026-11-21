#!/usr/bin/env python3
"""Step 5c — conversion or expansion?

For each disordered region placed in a homologous segment (or as an
insertion/extension), compare the length of the segment in the protein that
carries the region (A) with the ortholog's segment (B):
  internal  A span = residues between A's core anchors; B span = b_gap
  terminal  A tail beyond its last/first core anchor; B tail = b_gap
'gain / region length' near 0 means the tract replaced existing residues
(compositional conversion); near 1 means the tract is net added sequence
(expansion/insertion within the segment).
"""
import argparse
import csv
import json
import os

import numpy as np

PLAS = {"pfalciparum", "pvivax", "pknowlesi", "pberghei"}
LOCS = ("internal_variable", "internal_insertion", "terminal_variable", "terminal_extension")


def rd(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def spans(r, mem):
    LA = mem[(r["a_species"], r["a_gene"])]
    LB = mem[(r["b_species"], r["b_gene"])]
    if r["location"].startswith("internal"):
        a = int(r["a_anchor_right"]) - int(r["a_anchor_left"]) - 1
    elif r["a_anchor_right"] and not r["a_anchor_left"]:
        a = int(r["a_anchor_right"]) - 1
    else:
        a = LA - int(r["a_anchor_left"])
    return a, int(r["b_gap"]), LA, LB


def summary(sel):
    d = np.array([a - b for _, a, b, _, _ in sel])
    rl = np.array([int(r["region_len"]) for r, *_ in sel])
    rest = np.array([(LA - a) - (LB - b) for _, a, b, LA, LB in sel])
    f = d / rl
    return {"n": len(sel), "median_span_gain_aa": float(np.median(d)),
            "gain_over_region_len_median": round(float(np.median(f)), 3),
            "gain_over_region_len_iqr": [round(float(np.percentile(f, q)), 3) for q in (25, 75)],
            "pct_A_longer": round(100 * float(np.mean(d > 0)), 1),
            "rest_of_protein_gain_median_aa": float(np.median(rest))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--placement", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rp = rd(os.path.join(a.placement, "region_partner.tsv"))
    mem = {(r["species"], r["gene"]): int(r["model_len"])
           for r in rd(os.path.join(a.placement, "members.tsv"))}
    rows = []
    for r in rp:
        if (r["a_species"] in PLAS and r["b_species"] in PLAS and r["location"] in LOCS
                and float(r["region_plddt"]) < 70):
            s = spans(r, mem)
            if s[0] > 0:
                rows.append((r,) + s)
    pf = [x for x in rows if x[0]["a_species"] == "pfalciparum"]
    nonpf = [x for x in rows if x[0]["a_species"] != "pfalciparum"]
    res = {
        "Pf_Nrich": summary([x for x in pf if x[0]["lcr_class"] == "N-rich"]),
        "Pf_Nrich_partner_not_Nrich": summary([x for x in pf if x[0]["lcr_class"] == "N-rich"
                                               and x[0]["b_seg_class"] != "N-rich"]),
        "Pf_Nrich_partner_Nrich": summary([x for x in pf if x[0]["lcr_class"] == "N-rich"
                                           and x[0]["b_seg_class"] == "N-rich"]),
        "nonPf_Nrich": summary([x for x in nonpf if x[0]["lcr_class"] == "N-rich"]),
        "nonPf_other_classes": summary([x for x in nonpf if x[0]["lcr_class"] != "N-rich"]),
    }
    for cl in ("GS-rich", "charged", "other"):
        res[f"Pf_{cl}"] = summary([x for x in pf if x[0]["lcr_class"] == cl])
    for sp in ("pvivax", "pknowlesi", "pberghei"):
        res[f"Pf_Nrich_vs_{sp}"] = summary([x for x in pf if x[0]["lcr_class"] == "N-rich"
                                            and x[0]["b_species"] == sp])
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "segment_length_check.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
