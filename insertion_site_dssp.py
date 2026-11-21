#!/usr/bin/env python3
"""Item 4a — insertion-site geometry with DSSP instead of CA-geometry SSE.

Same rows and null as insertion_site_null.py (Plasmodium pairs, internal_insertion,
disordered region), restricted to sites with >= 1 B residue between the anchors so the
anchors can be recovered from region_partner.tsv (b1 = b_seg_start - 2, b2 = b_seg_end,
0-based). Both assignments (CA-geometry and DSSP) are computed on the same rows and the
same random draws, plus per-residue agreement between the two methods.
"""
import argparse, csv, json, math, os, random, sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import structlib as S          # noqa: E402
import dssp_lite as D          # noqa: E402
from insertion_site_null import site_sse, PLAS   # noqa: E402


def z(obs, null, used):
    n_null = sum(null.values())
    out = {}
    for k in ('loop', 'element', 'boundary'):
        p = null.get(k, 0) / n_null
        o = obs.get(k, 0)
        sd = math.sqrt(used * p * (1 - p)) or 1
        out[k] = dict(obs=round(o / used, 3), null=round(p, 3), z=round((o - used * p) / sd, 2))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--placement', required=True)
    ap.add_argument('--status', required=True)
    ap.add_argument('--structures', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--draws', type=int, default=200)
    a = ap.parse_args()
    rng = random.Random(1)
    path = {}
    with open(a.status) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['file']:
                path[(r['species'], r['gene'])] = os.path.join(a.structures, r['file'])
    with open(os.path.join(a.placement, 'region_partner.tsv')) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter='\t')
                if r['location'] == 'internal_insertion' and r['a_species'] in PLAS
                and r['b_species'] in PLAS and float(r['region_plddt']) < 70 and r['b_gap'] != '']
    n_all = len(rows)
    rows = [r for r in rows if r['b_seg_start']]
    cache = {}
    obs = {'ca': Counter(), 'dssp': Counter()}
    null = {'ca': Counter(), 'dssp': Counter()}
    agree = Counter()
    used = 0
    for r in rows:
        key = (r['b_species'], r['b_gene'])
        if key not in cache:
            m = S.read_model(path[key])
            X, _ = D.read_backbone(path[key])
            if len(X) != len(m['seq']):
                cache[key] = None
                continue
            m['ca_sse'] = S.ca_sse(m['ca'])
            m['dssp'] = D.assign(X)
            for i, (c1, c2) in enumerate(zip(m['ca_sse'], m['dssp'])):
                if m['plddt'][i] >= 70:
                    agree[(c1, c2)] += 1
            cache[key] = m
        m = cache[key]
        if m is None:
            continue
        b1, b2 = int(r['b_seg_start']) - 2, int(r['b_seg_end'])
        g = b2 - b1 - 1
        L = len(m['seq'])
        ok = [j for j in range(L - g - 1) if m['plddt'][j] >= 70 and m['plddt'][j + g + 1] >= 70]
        if not ok or b2 >= L:
            continue
        used += 1
        obs['ca'][site_sse(m['ca_sse'], b1, b2)] += 1
        obs['dssp'][site_sse(m['dssp'], b1, b2)] += 1
        for _ in range(a.draws):
            j = rng.choice(ok)
            null['ca'][site_sse(m['ca_sse'], j, j + g + 1)] += 1
            null['dssp'][site_sse(m['dssp'], j, j + g + 1)] += 1
    tot = sum(agree.values())
    res = dict(n_rows_all=n_all, n_with_recoverable_anchors=len(rows), n_used=used,
               ca_geometry=z(obs['ca'], null['ca'], used), dssp=z(obs['dssp'], null['dssp'], used),
               per_residue_agreement_plddt70=round(sum(v for (x, y), v in agree.items() if x == y) / tot, 3),
               confusion_ca_vs_dssp={f'{x}->{y}': v for (x, y), v in sorted(agree.items())},
               n_models=sum(1 for v in cache.values() if v))
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    json.dump(res, open(a.out, 'w'), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
