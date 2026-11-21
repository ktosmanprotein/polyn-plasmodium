#!/usr/bin/env python3
"""Item 5 — is AAT purity decay in Pf poly-N tracts consistent with neutral mutation?

Model: an N codon's third position flips AAT <-> AAC as a two-state stationary process with
equilibrium AAT fraction p_eq and total flip rate k. A tract born pure (AAT fraction p0 ~ 1)
t ago has expected purity  p(t) = p_eq + (p0 - p_eq) exp(-k t).

Calibration from the genome (no tract data used):
  p_eq  = AAT / (AAT + AAC) at Pf N codons outside low-complexity regions (unmasked in the RBH
          alignments).
  k t   per species pair from N<->N aligned codons outside LCRs:
          D = fraction differing at codon position 3;  k*(2t) = -ln(1 - D / (2 p_eq (1 - p_eq)))
          node depth (k t) = half of that.
Prediction for an age class whose origin lies between node depth a and b (uniform): mean
p(t) over [a, b]. Compared with the observed mean AAT fraction per class (step 6 dating).
A purity that stays above prediction = homogenisation (slippage copying, conversion) or
selection for AAT; at prediction = neutral decay; below = faster than neutral.
"""
import argparse, json, math, re, collections, statistics

PAIRS = [('ppraefalciparum', 'FP'), ('preichenowi', 'FPR'), ('pbillcollinsi', 'FPRB'),
         ('pblacklocki', 'CladeB'), ('padleri', 'Lav'), ('pgaboni', 'Lav'),
         ('pvivax', 'LavP'), ('pknowlesi', 'LavP'), ('pberghei', 'LavP'), ('pyoelii', 'LavP'),
         ('pgallinaceum', 'LavP'), ('prelictum', 'LavP')]
CLASS_BOUNDS = [('Pf only', None, 'FP'), ('Pf+Pprf', 'FP', 'FPR'), ('+Pr', 'FPR', 'FPRB'),
                ('+Pbilc', 'FPRB', 'CladeB'), ('clade B (+Pblac)', 'CladeB', 'Lav'),
                ('all Laverania (+clade A)', 'Lav', 'LavP')]


def norm(i):
    i = i.split()[0].lstrip('>')
    i = re.sub(r'-p\d+$', '', i)
    i = re.sub(r'-t\d+_\d+$', '', i)
    i = re.sub(r'_t\d+$', '', i)
    i = re.sub(r'\.\d+$', '', i)
    return i


def read_fa(path):
    d, k, buf = {}, None, []
    for line in open(path):
        if line.startswith('>'):
            if k:
                d[k] = ''.join(buf)
            k, buf = norm(line), []
        else:
            buf.append(line.strip())
    if k:
        d[k] = ''.join(buf)
    return d


def walk(m8, cdsA, cdsB, stats):
    for line in open(m8):
        f = line.rstrip('\n').split('\t')
        if len(f) < 13:
            continue
        qa, ta = norm(f[0]), norm(f[1])
        qs, ts = int(f[2]), int(f[4])
        qal, tal = f[11], f[12]
        A, B = cdsA.get(qa), cdsB.get(ta)
        if not A or not B:
            stats['missing_cds'] += 1
            continue
        i, j = qs - 1, ts - 1
        for a, b in zip(qal, tal):
            if a != '-' and b != '-':
                if a == 'N' and b == 'N':
                    ca, cb = A[3 * i:3 * i + 3].upper(), B[3 * j:3 * j + 3].upper()
                    if ca in ('AAT', 'AAC') and cb in ('AAT', 'AAC'):
                        stats['nn'] += 1
                        stats['diff'] += ca != cb
                        stats['pf_aat'] += ca == 'AAT'
                    else:
                        stats['codon_mismatch'] += 1
            if a != '-':
                i += 1
            if b != '-':
                j += 1


def mean_exp(a, b):
    if b == float('inf'):
        b = 50.0
    if b - a < 1e-9:
        return math.exp(-a)
    return (math.exp(-a) - math.exp(-b)) / (b - a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seqs', required=True)
    ap.add_argument('--rbh', required=True)
    ap.add_argument('--dating', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    pf = read_fa(f'{a.seqs}/pfalciparum_cds.fasta')
    pair_res = {}
    pf_aat = pf_n = 0
    for sp, node in PAIRS:
        cds = read_fa(f'{a.seqs}/{sp}_cds.fasta')
        st = collections.Counter()
        walk(f'{a.rbh}/rbh_pfalciparum__{sp}.m8', pf, cds, st)
        pair_res[sp] = dict(node=node, **st)
        pf_aat += st['pf_aat']
        pf_n += st['nn']
    # p_eq from Pf N codons aligned (outside LCRs) in the Laverania pairs, which cover most genes
    lav = [s for s, n in PAIRS if n in ('FP', 'FPR')]
    p_eq = sum(pair_res[s]['pf_aat'] for s in lav) / sum(pair_res[s]['nn'] for s in lav)
    H = 2 * p_eq * (1 - p_eq)
    depth = collections.defaultdict(list)
    for sp, r in pair_res.items():
        D = r['diff'] / r['nn'] if r['nn'] else float('nan')
        r['D'] = round(D, 4)
        if D < H * 0.98:
            kt2 = -math.log(1 - D / H)
            r['kt_node'] = round(kt2 / 2, 4)
        else:
            r['kt_node'] = float('inf')           # saturated
        depth[r['node']].append(r['kt_node'])
    node_kt = {n: (statistics.median(v) if all(math.isfinite(x) for x in v) else
                   (statistics.median([x for x in v if math.isfinite(x)]) if any(math.isfinite(x) for x in v) else float('inf')))
               for n, v in depth.items()}

    obs = collections.defaultdict(list)
    for line in open(a.dating):
        if line.startswith('gene'):
            hdr = line.rstrip('\n').split('\t')
            continue
        d = dict(zip(hdr, line.rstrip('\n').split('\t')))
        if d['aat_frac']:
            obs[d['origin']].append(float(d['aat_frac']))
    young = obs['Pf only']
    p0 = statistics.mean(young)
    rows = []
    for cls, lo, hi in CLASS_BOUNDS:
        a_ = 0.0 if lo is None else node_kt[lo]
        b_ = node_kt[hi]
        pred1 = p_eq + (1 - p_eq) * mean_exp(a_, b_)
        pred0 = p_eq + (p0 - p_eq) * mean_exp(a_, b_)
        v = obs.get(cls, [])
        m = statistics.mean(v)
        se = statistics.stdev(v) / math.sqrt(len(v))
        rows.append(dict(cls=cls, kt_from=round(a_, 4), kt_to=b_ if not math.isfinite(b_) else round(b_, 4),
                         n=len(v), observed=round(m, 3), se=round(se, 3),
                         predicted_born_pure=round(pred1, 3), predicted_born_at_young_mean=round(pred0, 3)))
    # rate needed to explain the observed old-class purity, if born pure, relative to neutral
    old = rows[-1]
    out = dict(p_eq=round(p_eq, 4), equilibrium_heterozygosity=round(H, 4), p0_young_mean=round(p0, 3),
               pairs=pair_res, node_depth_kt=node_kt, classes=rows)
    # per-class implied effective k t (if born pure): solve p = p_eq + (1-p_eq) e^{-x}
    for r in rows:
        x = (r['observed'] - p_eq) / (1 - p_eq)
        r['implied_kt_if_born_pure'] = round(-math.log(x), 4) if 0 < x < 1 else None
    json.dump(out, open(a.out + '.json', 'w'), indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != 'pairs'}, indent=1, default=str))
    for sp, r in pair_res.items():
        print(sp, r['node'], r['nn'], r['D'], r['kt_node'], 'mismatch', r.get('codon_mismatch', 0), 'missing', r.get('missing_cds', 0))


if __name__ == '__main__':
    main()
