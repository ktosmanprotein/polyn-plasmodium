#!/usr/bin/env python3
"""Step 6 — date P. falciparum poly-N tracts across Laverania and outgroups.

Inputs per species: primary-isoform proteins (headers = gene ids), fLPS
regions (summarise_lcr-style <key>_lcr.tsv) and MMseqs2 reciprocal-best-hit
alignments P. falciparum vs that species (qaln/taln).

For every P. falciparum N-rich fLPS region (the tract) and every species X:
  - the tract is placed in X's ortholog between flanking ANCHORS: the nearest
    aligned position on each side of the tract that starts/ends a run of
    >= 3 consecutive aligned (non-gap) columns. Residues inside the tract are
    never used as anchors, so poly-N is not placed by aligning N to N.
  - state:  present        X's segment between the anchors overlaps one of X's
                           own N-rich fLPS regions by >= 10 residues
            absent         X has a segment there but no N-rich region
            absent_noseg   X's segment is <= 5 residues (no counterpart)
            no_anchor      tract outside the alignment on one side (terminal);
                           X's tail beyond the last anchor is used as segment
            no_ortholog    no reciprocal best hit
  Terminal tracts use the tail and the same present/absent rule.

Dating (Dollo-style): a tract's origin is the most recent common ancestor of
P. falciparum and every species where it is present; species with no
ortholog are uninformative. Laverania topology after Otto et al. 2018:
((((Pf, Pprf), Pr), Pbilc), Pblac) in clade B; (Padl, Pgab) clade A.
Presence outside Laverania is reported by group without assuming where
Laverania attach (their position relative to avian vs mammalian parasites is
contested).
Homoplasy: absences inside the origin clade (requires extra losses or
independent gains) are counted per tract.
Codon purity: fraction of the tract's Asn codons that are AAT and the
longest pure AAT run, from the P. falciparum CDS.
"""
import argparse
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polylib as P  # noqa: E402

LAV_B = ["ppraefalciparum", "preichenowi", "pbillcollinsi", "pblacklocki"]
LAV_A = ["padleri", "pgaboni"]
OUT = {"primate": ["pvivax", "pknowlesi"], "rodent": ["pberghei", "pyoelii"],
       "avian": ["pgallinaceum", "prelictum"]}
SHORT = {"pfalciparum": "Pf", "ppraefalciparum": "Pprf", "preichenowi": "Pr",
         "pbillcollinsi": "Pbilc", "pblacklocki": "Pblac", "padleri": "Padl",
         "pgaboni": "Pgab", "pvivax": "Pv", "pknowlesi": "Pk", "pberghei": "Pb",
         "pyoelii": "Py", "pgallinaceum": "Pgal", "prelictum": "Prel"}
# clade B nested order: origin levels from youngest to oldest
LEVELS = [("Pf only", []),
          ("Pf+Pprf", ["ppraefalciparum"]),
          ("+Pr", ["preichenowi"]),
          ("+Pbilc", ["pbillcollinsi"]),
          ("clade B (+Pblac)", ["pblacklocki"]),
          ("all Laverania (+clade A)", LAV_A)]


def rd(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def aligned_map(h):
    """Pf position -> X position (0-based) for columns aligned in both."""
    m = {}
    qi, ti = int(h["qstart"]) - 1, int(h["tstart"]) - 1
    for a, b in zip(h["qaln"], h["taln"]):
        # masked (X) columns are low-complexity residues: never used as
        # alignment evidence, so anchors come only from ordinary sequence
        if a != "-" and b != "-" and a != "X" and b != "X":
            m[qi] = ti
        if a != "-":
            qi += 1
        if b != "-":
            ti += 1
    return m


def anchors(amap, s, e, run=3):
    """Nearest anchor left of s and right of e (0-based, exclusive e), each
    belonging to a run of >= `run` consecutive aligned Pf positions whose X
    partners are also consecutive."""
    def good(i, step):
        for k in range(run):
            j = i + step * k
            if j not in amap or amap[j] != amap[i] + step * k:
                return False
        return True
    left = next((i for i in range(s - 1, -1, -1) if i in amap and good(i, -1)), None)
    right = next((i for i in range(e, max(amap) + 1 if amap else e) if i in amap and good(i, 1)), None)
    return left, right


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True, help="folder with seqs/, lcr/, rbh/")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-overlap", type=int, default=10)
    ap.add_argument("--lenient", action="store_true",
                    help="count present_weak (segment >= 40%% N, below fLPS significance) as present")
    a = ap.parse_args()

    species = ["pfalciparum"] + LAV_B + LAV_A + sum(OUT.values(), [])
    prot, cds = {}, {}
    for sp in species:
        recs = [(P.parse_veupath_header(h), s.upper().rstrip("*"))
                for h, s in P.read_fasta(os.path.join(a.work, "seqs", f"{sp}_proteins.fasta"))]
        recs = [(m, s) for m, s in recs if not P.is_pseudo(m)]
        prim = P.primary_isoforms(recs)
        prot[sp] = {g: s for g, (_, s) in prim.items()}
        if sp == "pfalciparum":
            tx = {m["transcript"]: g for g, (m, _) in prim.items()}
            for h, s in P.read_fasta(os.path.join(a.work, "seqs", "pfalciparum_cds.fasta")):
                m = P.parse_veupath_header(h)
                g = tx.get(m["transcript"]) or tx.get(m["id"])
                if g:
                    cds[g] = s.upper()

    nrich = defaultdict(list)
    for sp in species:
        for r in rd(os.path.join(a.work, "lcr", f"{sp}_lcr.tsv")):
            if r["lcr_class"] == "N-rich":
                nrich[(sp, r["gene"])].append((int(r["start"]) - 1, int(r["end"])))

    rbh = {}
    for sp in species[1:]:
        p = os.path.join(a.work, "rbh", f"rbh_pfalciparum__{sp}.m8")
        rbh[sp] = {}
        if not os.path.exists(p):
            print(f"WARNING: {p} missing; {sp} treated as having no orthologs")
            continue
        cols = "query,target,qstart,qend,tstart,tend,qlen,tlen,pident,evalue,bits,qaln,taln".split(",")
        with open(p) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) == len(cols):
                    d = dict(zip(cols, f))
                    rbh[sp][d["query"]] = d

    rows = []
    for (sp, g), regs in sorted(nrich.items()):
        if sp != "pfalciparum" or g not in prot["pfalciparum"]:
            continue
        seq = prot["pfalciparum"][g]
        c = cds.get(g, "")
        for s, e in regs:
            tract = seq[s:e]
            codons = [c[3 * i:3 * i + 3] for i in range(s, e)] if len(c) >= 3 * e else []
            ncod = [x for x in codons if x in ("AAT", "AAC")]
            pure = longest = 0
            for x in codons:
                pure = pure + 1 if x == "AAT" else 0
                longest = max(longest, pure)
            row = {"gene": g, "start": s + 1, "end": e, "length": e - s,
                   "n_frac": round(tract.count("N") / (e - s), 3),
                   "aat_frac": round(sum(x == "AAT" for x in ncod) / len(ncod), 3) if ncod else "",
                   "longest_pure_aat": longest if codons else ""}
            for x in species[1:]:
                h = rbh[x].get(g)
                if not h:
                    row[SHORT[x]] = "no_ortholog"
                    continue
                amap = aligned_map(h)
                xg = h["target"]
                L = len(prot[x].get(xg, ""))
                if not amap or not L:
                    row[SHORT[x]] = "no_ortholog"
                    continue
                left, right = anchors(amap, s, e)
                if left is not None and right is not None:
                    xs, xe = amap[left] + 1, amap[right]
                elif left is None and right is not None:
                    xs, xe = 0, amap[right]
                elif right is None and left is not None:
                    xs, xe = amap[left] + 1, L
                else:
                    row[SHORT[x]] = "no_anchor"
                    continue
                ov = 0
                for rs, re_ in nrich.get((x, xg), []):
                    ov = max(ov, min(xe, re_) - max(xs, rs))
                xseg = prot[x][xg][xs:xe]
                if ov >= a.min_overlap:
                    st = "present"
                elif len(xseg) >= 8 and xseg.count("N") / len(xseg) >= 0.4:
                    st = "present_weak"     # N-rich but below fLPS significance
                elif xe - xs <= 5:
                    st = "absent_noseg"
                else:
                    st = "absent"
                row[SHORT[x]] = st
                row[SHORT[x] + "_seglen"] = max(xe - xs, 0)
            rows.append(row)

    present_states = {"present", "present_weak"} if a.lenient else {"present"}

    def is_p(r, sp):
        return r.get(SHORT[sp]) in present_states

    def is_a(r, sp):
        st = r.get(SHORT[sp], "")
        return st.startswith("absent") or (st == "present_weak" and not a.lenient)

    summary = Counter()
    homoplasy = Counter()
    by_age = defaultdict(list)
    out_presence = Counter()
    for r in rows:
        lav = LAV_B + LAV_A
        if not any(is_p(r, x) or is_a(r, x) for x in lav):
            r["origin"] = "uninformative"
        else:
            origin, members = "Pf only", []
            for name, add in LEVELS:
                members = members + add
                if any(is_p(r, x) for x in add):
                    origin = name
            r["origin"] = origin
            # absences inside the origin clade
            inside = []
            for name, add in LEVELS:
                inside += add
                if name == origin:
                    break
            r["internal_absences"] = sum(is_a(r, x) for x in inside)
            homoplasy[(origin, r["internal_absences"] > 0)] += 1
        for grp, sps in OUT.items():
            if any(is_p(r, x) for x in sps):
                out_presence[(r["origin"], grp)] += 1
        r["outside_laverania"] = ",".join(g for g, sps in OUT.items() if any(is_p(r, x) for x in sps))
        summary[r["origin"]] += 1
        by_age[r["origin"]].append(r)

    order = [n for n, _ in LEVELS] + ["uninformative"]
    res = {"n_tracts": len(rows), "origin_counts": {k: summary[k] for k in order},
           "clade_consistent_fraction": {}, "outside_laverania_by_origin": {},
           "tract_properties_by_origin": {}}
    for k in order[:-1]:
        n = summary[k]
        if n:
            res["clade_consistent_fraction"][k] = round(homoplasy[(k, False)] / n, 3)
            res["outside_laverania_by_origin"][k] = {g: out_presence[(k, g)] for g in OUT}
            v = by_age[k]
            def med(key):
                xs = [float(x[key]) for x in v if x[key] != ""]
                return round(float(np.median(xs)), 3) if xs else None
            res["tract_properties_by_origin"][k] = {
                "n": n, "median_length": med("length"), "median_n_frac": med("n_frac"),
                "median_aat_frac": med("aat_frac"), "median_longest_pure_aat": med("longest_pure_aat")}
    states = {SHORT[x]: dict(Counter(r.get(SHORT[x]) for r in rows)) for x in species[1:]}
    res["state_counts_by_species"] = states

    os.makedirs(a.out, exist_ok=True)
    fields = ["gene", "start", "end", "length", "n_frac", "aat_frac", "longest_pure_aat", "origin",
              "internal_absences", "outside_laverania"] + \
             [SHORT[x] for x in species[1:]] + [SHORT[x] + "_seglen" for x in species[1:]]
    tag = "_lenient" if a.lenient else ""
    P.write_tsv(os.path.join(a.out, f"pf_tract_dating{tag}.tsv"), rows, fields)
    with open(os.path.join(a.out, f"pf_tract_dating{tag}.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
