#!/usr/bin/env python3
"""Item 4b — second low-complexity caller (SEG) vs fLPS.

SEG after Wootton & Federhen (1993/1996), standard protein parameters W = 12, K1 = 2.2,
K2 = 2.5 bits: windows with complexity <= K1 trigger; triggered windows are merged and
extended while flanking windows stay <= K2. (The final optimal-subsequence trimming of SEG is
omitted, so segments can be a few residues longer than NCBI seg output.)
SEG segments >= 15 aa are classed like fLPS regions by their dominant residues
(N-rich if N >= 30% and N is the most common residue; Q-rich; GS-rich; charged (KEDR >= 50%);
other).

Outputs per species: SEG segment counts by class, % proteins with an N-rich SEG segment
(compared with fLPS), and for Pf: overlap of fLPS N-rich tracts with SEG N-rich segments
both ways.
"""
import argparse, csv, json, math, os, sys, collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polylib as P   # noqa: E402

W, K1, K2 = 12, 2.2, 2.5
SPECIES = ['pfalciparum', 'ppraefalciparum', 'preichenowi', 'pbillcollinsi', 'pblacklocki',
           'padleri', 'pgaboni', 'pvivax', 'pknowlesi', 'pberghei', 'pyoelii', 'pgallinaceum',
           'prelictum']


def entropy(counts, w):
    h = 0.0
    for c in counts.values():
        if c:
            p = c / w
            h -= p * math.log2(p)
    return h


def seg(seq):
    n = len(seq)
    if n < W:
        return []
    comp = []
    c = collections.Counter(seq[:W])
    comp.append(entropy(c, W))
    for i in range(1, n - W + 1):
        c[seq[i - 1]] -= 1
        c[seq[i + W - 1]] += 1
        comp.append(entropy(c, W))
    segs = []
    i = 0
    m = len(comp)
    while i < m:
        if comp[i] <= K1:
            s = i
            e = i
            while e + 1 < m and comp[e + 1] <= K2:
                e += 1
            while s - 1 >= 0 and comp[s - 1] <= K2:
                s -= 1
            segs.append((s, e + W))       # residues [s, e+W)
            i = e + 1
        else:
            i += 1
    merged = []
    for s, e in segs:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(e, merged[-1][1]))
        else:
            merged.append((s, e))
    return merged


def classify(sub):
    n = len(sub)
    c = collections.Counter(sub)
    top, _ = c.most_common(1)[0]
    if c['N'] / n >= 0.3 and top == 'N':
        return 'N-rich'
    if c['Q'] / n >= 0.3 and top == 'Q':
        return 'Q-rich'
    if (c['G'] + c['S']) / n >= 0.4 and top in 'GS':
        return 'GS-rich'
    if sum(c[x] for x in 'KEDR') / n >= 0.5:
        return 'charged'
    return 'other'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seqs', required=True)
    ap.add_argument('--lcr', required=True, nargs='+', help='folders with <species>_lcr.tsv')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    res = {}
    pf_seg_n, pf_flps_n = [], []
    for sp in SPECIES:
        recs = [(P.parse_veupath_header(h), s.upper().rstrip('*'))
                for h, s in P.read_fasta(os.path.join(a.seqs, f'{sp}_proteins.fasta'))]
        recs = [(m, s) for m, s in recs if not P.is_pseudo(m)]
        prim = P.primary_isoforms(recs)
        cls_counts = collections.Counter()
        prot_n = set()
        segN = collections.defaultdict(list)
        for g, (_, s) in prim.items():
            for st, en in seg(s):
                if en - st < 15:
                    continue
                k = classify(s[st:en])
                cls_counts[k] += 1
                if k == 'N-rich':
                    prot_n.add(g)
                    segN[g].append((st, en))
        flN = collections.defaultdict(list)
        for d in a.lcr:
            p = os.path.join(d, f'{sp}_lcr.tsv')
            if os.path.exists(p):
                with open(p) as fh:
                    for r in csv.DictReader(fh, delimiter='\t'):
                        if r['lcr_class'] == 'N-rich':
                            flN[r['gene']].append((int(r['start']) - 1, int(r['end'])))
                break
        nprot = len(prim)
        res[sp] = dict(n_proteins=nprot, seg_segments_by_class=dict(cls_counts),
                       pct_proteins_Nrich_seg=round(100 * len(prot_n) / nprot, 1),
                       pct_proteins_Nrich_flps=round(100 * len(flN) / nprot, 1) if flN else None)

        def ov(A, B):
            hit = 0
            for g, iv in A.items():
                for s, e in iv:
                    if any(min(e, e2) - max(s, s2) >= 5 for s2, e2 in B.get(g, [])):
                        hit += 1
            tot = sum(len(v) for v in A.values())
            return hit, tot
        if flN:
            h1, t1 = ov(flN, segN)
            h2, t2 = ov(segN, flN)
            res[sp]['flps_Nrich_found_by_seg'] = [h1, t1, round(h1 / t1, 3) if t1 else None]
            res[sp]['seg_Nrich_found_by_flps'] = [h2, t2, round(h2 / t2, 3) if t2 else None]
            if sp == 'pfalciparum':
                miss = [(g, s, e) for g, iv in segN.items() for s, e in iv
                        if not any(min(e, e2) - max(s, s2) >= 5 for s2, e2 in flN.get(g, []))]
                res[sp]['seg_only_Nrich_median_len'] = sorted(e - s for _, s, e in miss)[len(miss) // 2] if miss else None
    json.dump(res, open(a.out, 'w'), indent=1)
    for sp, r in res.items():
        print(sp, r['pct_proteins_Nrich_seg'], r['pct_proteins_Nrich_flps'], r.get('flps_Nrich_found_by_seg'), r.get('seg_Nrich_found_by_flps'))


if __name__ == '__main__':
    main()
