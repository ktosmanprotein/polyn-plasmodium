#!/usr/bin/env python3
"""Check 2 — DNA mismatch-repair (MMR) and replication genes across the 13 species.

For each Pf gene in the set: reciprocal-best-hit ortholog in each species (step-6 RBH
alignments), Pf length covered by the alignment, identity, and the identity's percentile among
all RBH pairs for that species (genome-wide background). A gene lost, truncated, or unusually
divergent on the Laverania side would be a candidate cause of the Laverania tract burst.
Control pairs that do not involve Laverania (Pv-Pb, Pk-Pb, Pv-Pk; step-8 rbh_pilot) show whether
divergence is Laverania-specific. Product names in each species' protein headers are searched as
a second route, so RBH failures are not mistaken for gene loss. Key catalytic motifs are checked
in every ortholog: MutS Walker A (G[PSA]N..GK[ST]) for MSH proteins; MutL/GHKL motif I
(G[FY]RGE[AG]L) for MLH/PMS proteins.
"""
import argparse, csv, json, os, re, sys, collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polylib as P   # noqa: E402

GENES = {
    'PF3D7_0706700': 'MSH2-1', 'PF3D7_1427500': 'MSH2-2', 'PF3D7_0505500': 'MSH6',
    'PF3D7_1117800': 'MLH1', 'PF3D7_0726300': 'PMS1', 'PF3D7_1405400': 'MMR (MLH-like)',
    'PF3D7_1017000': 'POLD1', 'PF3D7_0308000': 'POLD2', 'PF3D7_0630300': 'POLE1',
    'PF3D7_1234300': 'POLE2', 'PF3D7_1361900': 'PCNA1', 'PF3D7_1226600': 'PCNA2',
}
MOTIF = {'MSH': re.compile(r'G[PSA]N..GK[ST]'), 'MLH': re.compile(r'G[FY]RGE[AG]L')}
SPECIES = ['ppraefalciparum', 'preichenowi', 'pbillcollinsi', 'pblacklocki', 'padleri', 'pgaboni',
           'pvivax', 'pknowlesi', 'pberghei', 'pyoelii', 'pgallinaceum', 'prelictum']
COL = "query,target,qstart,qend,tstart,tend,qlen,tlen,pident,evalue,bits,qaln,taln".split(",")


def motif_kind(name):
    return 'MSH' if name.startswith('MSH') else ('MLH' if name.startswith(('MLH', 'PMS', 'MMR')) else None)


def load_prot(seqs, sp):
    recs = [(P.parse_veupath_header(h), s.upper().rstrip('*'))
            for h, s in P.read_fasta(os.path.join(seqs, f'{sp}_proteins.fasta'))]
    recs = [(m, s) for m, s in recs if not P.is_pseudo(m)]
    prim = P.primary_isoforms(recs)
    return {g: s for g, (_, s) in prim.items()}, {g: m for g, (m, _) in prim.items()}


def rbh_table(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 11:
                rows.append(dict(zip(COL, f)))
    return rows


def pid(r):
    """Identity over aligned columns where neither residue is masked (X) or a gap; falls back
    to the reported pident when no alignment strings are present."""
    if r.get('qaln') and r.get('taln'):
        m = n = 0
        for x, y in zip(r['qaln'], r['taln']):
            if x in '-X' or y in '-X':
                continue
            n += 1
            m += x == y
        return 100 * m / n if n else 0.0
    v = float(r['pident'])
    return v * 100 if v <= 1 else v


def pct_rank(vals, x):
    return round(100 * sum(v <= x for v in vals) / len(vals), 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seqs', required=True)
    ap.add_argument('--rbh', required=True)
    ap.add_argument('--pilot', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    pf, pfmeta = load_prot(a.seqs, 'pfalciparum')
    res = {g: {'name': n, 'pf_len': len(pf.get(g, '')),
               'pf_motif': bool(MOTIF[motif_kind(n)].search(pf[g])) if motif_kind(n) else None,
               'species': {}} for g, n in GENES.items()}
    bg = {}
    for sp in SPECIES:
        prot, meta = load_prot(a.seqs, sp)
        rows = rbh_table(os.path.join(a.rbh, f'rbh_pfalciparum__{sp}.m8'))
        allp = [pid(r) for r in rows]
        bg[sp] = dict(n=len(allp), median_pident=round(sorted(allp)[len(allp) // 2], 1))
        by = {r['query']: r for r in rows}
        prod = collections.defaultdict(list)
        for g, m in meta.items():
            prod[(m.get('gene_product') or m.get('product') or '').lower()].append(g)
        for g, n in GENES.items():
            r = by.get(g)
            ent = {}
            if r:
                t = r['target']
                cov = (int(r['qend']) - int(r['qstart']) + 1) / max(1, int(r['qlen']))
                ent = dict(ortholog=t, len=len(prot.get(t, '')), pf_cov=round(cov, 2), pident=round(pid(r), 1),
                           pident_percentile=pct_rank(allp, pid(r)))
                k = motif_kind(n)
                if k and t in prot:
                    ent['motif'] = bool(MOTIF[k].search(prot[t]))
            else:
                key = (pfmeta[g].get('gene_product') or '').lower()
                ent = dict(ortholog=None, same_product_genes=prod.get(key, [])[:5])
            res[g]['species'][sp] = ent
    # control pairs not involving Laverania
    ctrl = {}
    for f in ['rbh_pvivax__pberghei.tsv', 'rbh_pknowlesi__pberghei.tsv', 'rbh_pvivax__pknowlesi.tsv']:
        p = os.path.join(a.pilot, f)
        if not os.path.exists(p):
            continue
        rows = rbh_table(p)
        allp = [pid(r) for r in rows]
        ctrl[f] = dict(n=len(allp), median_pident=round(sorted(allp)[len(allp) // 2], 1), genes={})
        # map Pf gene -> Pv/Pk ortholog via step-6 RBH, then look up control pair
        q = f.split('__')[0].replace('rbh_', '')
        pf2q = {r['query']: r['target'] for r in rbh_table(os.path.join(a.rbh, f'rbh_pfalciparum__{q}.m8'))}
        byq = {r['query']: r for r in rows}
        for g, n in GENES.items():
            o = pf2q.get(g)
            r = byq.get(o) if o else None
            ctrl[f]['genes'][n] = dict(pident=round(pid(r), 1), percentile=pct_rank(allp, pid(r))) if r else None
    out = dict(background=bg, genes=res, control_pairs=ctrl)
    json.dump(out, open(a.out + '.json', 'w'), indent=1)
    # compact table
    with open(a.out + '.tsv', 'w') as fh:
        fh.write('gene\tname\tpf_len\t' + '\t'.join(f'{s}_pid(pct)' for s in SPECIES) + '\n')
        for g, d in res.items():
            cells = []
            for s in SPECIES:
                e = d['species'][s]
                if e.get('ortholog'):
                    m = '' if e.get('motif') in (None, True) else ' noMotif'
                    cells.append(f"{e['pident']}({e['pident_percentile']}) cov{e['pf_cov']}{m}")
                else:
                    cells.append('NONE' + (f" [{','.join(e['same_product_genes'])}]" if e.get('same_product_genes') else ''))
            fh.write(f"{g}\t{d['name']}\t{d['pf_len']}\t" + '\t'.join(cells) + '\n')
    print(open(a.out + '.tsv').read())
    print(json.dumps(bg))
    print(json.dumps({k: v['genes'] for k, v in ctrl.items()}))


if __name__ == '__main__':
    main()
