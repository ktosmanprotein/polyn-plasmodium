#!/usr/bin/env python3
"""Population length variation of Pf poly-N tracts by tract age (Pf8, MalariaGEN).

Input
  tract_regions.json   2,202 Pf tracts {id, chrom, start, end, origin, len, aat, pure}
  pf8_tract_variants.tsv  site-level Pf8 records near each tract, extracted from
                       Pf8.zarr (variants/POS, CHROM, FILTER_PASS, AN, AC, REF, ALT)
                       columns: tract, pos, pass, AN, REF, ALTs(,), ACs(,)

Per tract (PASS sites overlapping [start-1, end] only):
  n_len        number of length-changing alleles with AC >= 1
  n_len_1pct   number with frequency >= 1%
  len_het      1 - sum(f_i^2) over length classes (ref + each indel allele), the
               chance two random genomes differ in tract length at the commonest site
               summed over sites (approximation: sites treated independently)
  max_len_f    highest single indel-allele frequency
  snp_per_bp   PASS SNV alleles (AC >= 1) per bp, as a mutability / callability control
  an_frac      median AN / max AN over sites (callability)
Then compares age classes with length/purity control (logistic and linear models).
"""
import argparse, json, math, collections
import numpy as np

AGE = ['Pf only', 'Pf+Pprf', '+Pr', '+Pbilc', 'clade B (+Pblac)', 'all Laverania (+clade A)']
SHORT = {'Pf only': 'Pf', 'Pf+Pprf': 'Pf+Pprf', '+Pr': '+Pr', '+Pbilc': '+Pbilc',
         'clade B (+Pblac)': 'cladeB', 'all Laverania (+clade A)': 'allLav'}


def load_variants(path):
    by = collections.defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith('tract'):
                continue
            t, pos, ps, an, ref, alts, acs = line.rstrip('\n').split('\t')
            by[int(t)].append((int(pos), ps == '1', int(an), ref,
                               alts.split(',') if alts else [],
                               [int(x) for x in acs.split(',')] if acs else []))
    return by


def tract_metrics(tr, recs, max_an):
    s, e = tr['start'], tr['end']
    n_len = n_len1 = n_fs = n_inf = 0
    het = 0.0
    maxf = 0.0
    snp = 0
    ans = []
    nsite = 0
    for pos, ps, an, ref, alts, acs in recs:
        if not ps or an <= 0:
            continue
        if pos > e or pos + len(ref) - 1 < s - 1:
            continue
        nsite += 1
        ans.append(an)
        fsum2 = 0.0
        fref = 1.0
        for a, c in zip(alts, acs):
            if a == '*' or c <= 0:
                continue
            f = c / an
            if len(a) != len(ref):
                d = abs(len(a) - len(ref))
                n_len += 1
                if d % 3 == 0:
                    n_inf += 1
                else:
                    n_fs += 1
                if f >= 0.01:
                    n_len1 += 1
                maxf = max(maxf, f)
                fsum2 += f * f
                fref -= f
            elif len(a) == 1 and len(ref) == 1:
                snp += 1
        fref = max(fref, 0.0)
        het += 1 - (fsum2 + fref * fref)
    span = e - s + 1
    return dict(id=tr['id'], origin=tr['origin'], len=tr['len'], aat=float(tr['aat']),
                pure=int(tr['pure']), span=span, nsite=nsite, n_len=n_len, n_len_1pct=n_len1,
                n_inframe=n_inf, n_frameshift=n_fs, len_het=het, max_len_f=maxf,
                snp_per_bp=snp / span, an_frac=(float(np.median(ans)) / max_an) if ans else float('nan'))


def logit(X, y, iters=50):
    X = np.column_stack([np.ones(len(y)), X])
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-X @ b))
        W = p * (1 - p)
        H = X.T @ (X * W[:, None]) + 1e-9 * np.eye(X.shape[1])
        g = X.T @ (y - p)
        step = np.linalg.solve(H, g)
        b += step
        if np.abs(step).max() < 1e-8:
            break
    cov = np.linalg.inv(H)
    se = np.sqrt(np.diag(cov))
    return b, se


def ols(X, y):
    X = np.column_stack([np.ones(len(y)), X])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    s2 = res @ res / (len(y) - X.shape[1])
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    return b, se


def pnorm(z):
    return math.erfc(abs(z) / math.sqrt(2))


def mwu(a, b):
    from scipy.stats import mannwhitneyu
    return float(mannwhitneyu(a, b).pvalue)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--regions', required=True)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    regs = json.load(open(a.regions))
    var = load_variants(a.variants)
    max_an = max(r[2] for v in var.values() for r in v)
    rows = [tract_metrics(t, var.get(t['id'], []), max_an) for t in regs]
    rows = [r for r in rows if r['origin'] in AGE]

    # spliced tracts (span >> 3*len) inflate span; keep tracts whose span <= 3*len + 30
    ok = [r for r in rows if r['span'] <= 3 * r['len'] + 30]

    summary = {}
    for g in AGE:
        rr = [r for r in ok if r['origin'] == g]
        if not rr:
            continue
        summary[SHORT[g]] = dict(
            n=len(rr), median_len=float(np.median([r['len'] for r in rr])),
            median_aat=float(np.median([r['aat'] for r in rr])),
            frac_any_len=float(np.mean([r['n_len'] > 0 for r in rr])),
            frac_len_1pct=float(np.mean([r['n_len_1pct'] > 0 for r in rr])),
            mean_n_len=float(np.mean([r['n_len'] for r in rr])),
            mean_len_het=float(np.mean([r['len_het'] for r in rr])),
            median_len_het=float(np.median([r['len_het'] for r in rr])),
            frac_inframe=(sum(r['n_inframe'] for r in rr) / max(1, sum(r['n_len'] for r in rr))),
            mean_snp_per_bp=float(np.mean([r['snp_per_bp'] for r in rr])),
            median_an_frac=float(np.nanmedian([r['an_frac'] for r in rr])))

    # young (Pf, Pf+Pprf) vs old (allLav)
    young = [r for r in ok if r['origin'] in AGE[:2]]
    mid = [r for r in ok if r['origin'] in AGE[2:5]]
    old = [r for r in ok if r['origin'] == AGE[5]]
    tests = {}
    for m in ['len_het', 'n_len', 'snp_per_bp']:
        tests[m + '_young_vs_old_p'] = mwu([r[m] for r in young], [r[m] for r in old])

    # models: response ~ age_rank + log(len) + aat
    rank = {g: i for i, g in enumerate(AGE)}   # 0 = youngest
    X = np.array([[rank[r['origin']], math.log(r['len']), r['aat']] for r in ok])
    y1 = np.array([r['n_len_1pct'] > 0 for r in ok], float)
    b, se = logit(X, y1)
    tests['logit_any1pct'] = {k: dict(coef=float(b[i]), se=float(se[i]), p=pnorm(b[i] / se[i]))
                             for i, k in enumerate(['const', 'age_rank', 'log_len', 'aat'])}
    y2 = np.log1p(np.array([r['len_het'] for r in ok]) * 100)
    b, se = ols(X, y2)
    tests['ols_log_len_het'] = {k: dict(coef=float(b[i]), se=float(se[i]), p=pnorm(b[i] / se[i]))
                               for i, k in enumerate(['const', 'age_rank', 'log_len', 'aat'])}
    # length-matched young vs old (bin by len)
    bins = [0, 12, 16, 22, 32, 1000]
    lm = []
    for lo, hi in zip(bins, bins[1:]):
        yy = [r for r in young if lo <= r['len'] < hi]
        oo = [r for r in old if lo <= r['len'] < hi]
        if yy and oo:
            lm.append(dict(bin=f'{lo}-{hi}', n_young=len(yy), n_old=len(oo),
                           young_any1pct=float(np.mean([r['n_len_1pct'] > 0 for r in yy])),
                           old_any1pct=float(np.mean([r['n_len_1pct'] > 0 for r in oo])),
                           young_het=float(np.mean([r['len_het'] for r in yy])),
                           old_het=float(np.mean([r['len_het'] for r in oo]))))
    tests['length_matched'] = lm
    # purity: within all tracts, does AAT fraction / pure run predict length polymorphism?
    from scipy.stats import spearmanr
    tests['spearman_aat_vs_het'] = float(spearmanr([r['aat'] for r in ok], [r['len_het'] for r in ok]).correlation)
    tests['spearman_pure_vs_het'] = float(spearmanr([r['pure'] for r in ok], [r['len_het'] for r in ok]).correlation)
    tests['spearman_len_vs_het'] = float(spearmanr([r['len'] for r in ok], [r['len_het'] for r in ok]).correlation)

    out = dict(n_tracts=len(rows), n_unspliced=len(ok), max_AN=max_an, by_age=summary, tests=tests)
    json.dump(out, open(a.out + '.json', 'w'), indent=1)
    with open(a.out + '.tsv', 'w') as f:
        keys = list(rows[0].keys())
        f.write('\t'.join(keys) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[k]) for k in keys) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
