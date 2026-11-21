#!/usr/bin/env python3
"""Item 5b — direct within-tract test: codon change inside shared tracts vs genome background.

For Pf tracts that are also present in a close relative (Pprf, Pr), the ortholog's segment is
found between the same flanking anchors used in step 6 (never by aligning N to N). Then:
  * AAT fraction of N codons in the Pf tract vs in the ortholog segment: mean |difference|,
    compared with a point-mutation null that flips each Pf N codon's third position with the
    genome-wide N-codon divergence D for that pair (codon_decay.json), 1,000 replicates.
  * Equal-length tracts (same number of codons, no length change): codon-by-codon N<->N
    third-position divergence inside tracts vs the background D.
"""
import argparse, json, os, sys, random, statistics, collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from codon_decay import read_fa, norm          # noqa: E402
from date_tracts import aligned_map, anchors   # noqa: E402

COL = "query,target,qstart,qend,tstart,tend,qlen,tlen,pident,evalue,bits,qaln,taln".split(",")
SP = {'ppraefalciparum': 'Pprf', 'preichenowi': 'Pr'}


def codons(s):
    return [s[i:i + 3] for i in range(0, len(s) - len(s) % 3, 3)]


def aat_frac(cs):
    n = [c for c in cs if c in ('AAT', 'AAC')]
    return (sum(c == 'AAT' for c in n) / len(n), len(n)) if n else (None, 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seqs', required=True)
    ap.add_argument('--rbh', required=True)
    ap.add_argument('--dating', required=True)
    ap.add_argument('--decay', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    D = {sp: r['D'] for sp, r in json.load(open(a.decay))['pairs'].items()}
    pf = read_fa(f'{a.seqs}/pfalciparum_cds.fasta')
    rows = []
    with open(a.dating) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for line in f:
            rows.append(dict(zip(hdr, line.rstrip('\n').split('\t'))))
    rng = random.Random(1)
    out = {}
    for sp, short in SP.items():
        cds = read_fa(f'{a.seqs}/{sp}_cds.fasta')
        hits = {}
        for line in open(f'{a.rbh}/rbh_pfalciparum__{sp}.m8'):
            h = dict(zip(COL, line.rstrip('\n').split('\t')))
            hits[norm(h['query'])] = h
        res = []
        for r in rows:
            if r[short] != 'present':
                continue
            g = r['gene']
            h = hits.get(g)
            if not h or g not in pf or norm(h['target']) not in cds:
                continue
            s, e = int(r['start']) - 1, int(r['end'])
            amap = aligned_map(h)
            if not amap:
                continue
            L, R = anchors(amap, s, e)
            if L is None or R is None:
                continue
            xs, xe = amap[L] + 1, amap[R]
            pf_seg = codons(pf[g][3 * (L + 1):3 * R].upper())
            ox = codons(cds[norm(h['target'])][3 * xs:3 * xe].upper())
            fa, na = aat_frac(pf_seg)
            fb, nb = aat_frac(ox)
            if fa is None or fb is None or na < 5 or nb < 5:
                continue
            # point-mutation null for |delta|
            null = []
            for _ in range(200):
                sim = [('AAC' if c == 'AAT' else 'AAT') if (c in ('AAT', 'AAC') and rng.random() < D[sp]) else c
                       for c in pf_seg]
                fs, _ = aat_frac(sim)
                null.append(abs(fs - fa))
            same_len = len(pf_seg) == len(ox)
            nn = diff = 0
            if same_len:
                for c1, c2 in zip(pf_seg, ox):
                    if c1 in ('AAT', 'AAC') and c2 in ('AAT', 'AAC'):
                        nn += 1
                        diff += c1 != c2
            res.append(dict(gene=g, origin=r['origin'], pf_aat=fa, x_aat=fb, n_pf=na, n_x=nb,
                            delta=abs(fa - fb), null_delta=statistics.mean(null), same_len=same_len,
                            nn=nn, diff=diff))
        sl = [x for x in res if x['same_len']]
        nn = sum(x['nn'] for x in sl)
        df = sum(x['diff'] for x in sl)
        out[short] = dict(
            n_tracts=len(res), background_D=D[sp],
            mean_abs_delta_aat=round(statistics.mean(x['delta'] for x in res), 4),
            null_mean_abs_delta=round(statistics.mean(x['null_delta'] for x in res), 4),
            frac_tracts_delta_gt_null95=None,
            n_same_length=len(sl), same_length_nn=nn, same_length_diff=df,
            within_tract_D=round(df / nn, 4) if nn else None,
            ratio_within_vs_background=round((df / nn) / D[sp], 2) if nn else None,
            mean_pf_aat=round(statistics.mean(x['pf_aat'] for x in res), 3),
            mean_x_aat=round(statistics.mean(x['x_aat'] for x in res), 3),
            frac_length_changed=round(1 - len(sl) / len(res), 3))
        with open(f'{a.out}_{short}.tsv', 'w') as fh:
            ks = list(res[0].keys())
            fh.write('\t'.join(ks) + '\n')
            for x in res:
                fh.write('\t'.join(str(x[k]) for k in ks) + '\n')
    json.dump(out, open(a.out + '.json', 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
