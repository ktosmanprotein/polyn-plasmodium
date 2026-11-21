#!/usr/bin/env python3
"""Step 5b — are true insertion sites special?

For every disordered region placed as internal_insertion (ortholog B lacks
the segment), the site in B is the pair of core anchors b1, b2 flanking it.
Observed: B's secondary structure at the site (loop / element / boundary) and
the CA-CA distance b1-b2. Null: for the same B model and the same anchor
spacing (b_gap), sites drawn at random from positions where both anchors are
confidently modelled (pLDDT >= 70). The comparison asks whether insertions
land at different kinds of places than a random confidently modelled spot
in the same protein.
"""
import argparse
import csv
import json
import os
import random
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import structlib as S  # noqa: E402

PLAS = {"pfalciparum", "pvivax", "pknowlesi", "pberghei"}


def site_sse(sse, b1, b2):
    s1, s2 = sse[b1], sse[b2]
    if s1 == "C" and s2 == "C":
        return "loop"
    if s1 == s2 and s1 in "HE" and all(c == s1 for c in sse[b1:b2 + 1]):
        return "element"
    return "boundary"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--placement", required=True)
    ap.add_argument("--status", required=True)
    ap.add_argument("--structures", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)

    path = {}
    with open(a.status) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["file"]:
                path[(r["species"], r["gene"])] = os.path.join(a.structures, r["file"])
    with open(os.path.join(a.placement, "region_partner.tsv")) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["location"] == "internal_insertion" and r["a_species"] in PLAS
                and r["b_species"] in PLAS and float(r["region_plddt"]) < 70
                and r["b_gap"] != ""]

    cache = {}
    obs, null = Counter(), Counter()
    obs_cls = defaultdict(Counter)
    dist_obs, dist_null = [], []
    used = 0
    for r in rows:
        key = (r["b_species"], r["b_gene"])
        if key not in cache:
            m = S.read_model(path[key])
            m["sse"] = S.ca_sse(m["ca"])
            cache[key] = m
        m = cache[key]
        g = int(r["b_gap"])
        L = len(m["seq"])
        ok = [j for j in range(L - g - 1)
              if m["plddt"][j] >= 70 and m["plddt"][j + g + 1] >= 70]
        if not ok:
            continue
        used += 1
        obs[r["b_site_sse"]] += 1
        obs_cls[r["lcr_class"]][r["b_site_sse"]] += 1
        if r["b_anchor_dist"]:
            dist_obs.append(float(r["b_anchor_dist"]))
        for _ in range(a.draws):
            j = rng.choice(ok)
            null[site_sse(m["sse"], j, j + g + 1)] += 1
            dist_null.append(S.dist(m["ca"][j], m["ca"][j + g + 1]))

    n_null = sum(null.values())
    res = {
        "n_insertions": used,
        "observed": {k: [v, round(v / used, 3)] for k, v in obs.items()},
        "null_fraction": {k: round(v / n_null, 3) for k, v in null.items()},
        "observed_by_class": {k: dict(v) for k, v in obs_cls.items()},
        "anchor_dist_median_observed": sorted(dist_obs)[len(dist_obs) // 2] if dist_obs else None,
        "anchor_dist_median_null": sorted(dist_null)[len(dist_null) // 2] if dist_null else None,
    }
    # binomial-style z for each category vs null fraction
    import math
    z = {}
    for k in ("loop", "element", "boundary"):
        p = null.get(k, 0) / n_null
        o = obs.get(k, 0)
        sd = math.sqrt(used * p * (1 - p)) or 1
        z[k] = round((o - used * p) / sd, 2)
    res["z_vs_null"] = z
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "insertion_site_null.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
