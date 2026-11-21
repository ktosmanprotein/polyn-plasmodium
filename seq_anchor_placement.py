#!/usr/bin/env python3
"""Item 4c — sequence-anchor placement of LCRs, validated against structure-based placement,
then applied to proteins without AlphaFold models.

For every pilot-species pair with a reciprocal best hit (MMseqs2 on LCR-masked proteins;
masked X columns never count as alignment evidence) and every fLPS region of A:
  anchors = nearest aligned positions left/right of the region that start a run of >= 3
            consecutive aligned columns (same rule as step-6 dating)
  terminal          no anchor on one side
  insertion         B segment between anchors <= max(5, 0.25 * A span between anchors)
  opposite_segment  otherwise (B has a segment of comparable size at the site)
Validation: rows that structure placement (placement_v3/region_partner.tsv) called
internal_insertion / internal_variable / terminal_* are compared with the sequence call.
Application: regions in proteins where A or B has no AlphaFold model.
"""
import argparse, csv, json, os, sys, collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from date_tracts import aligned_map, anchors        # noqa: E402

PILOT = ['pfalciparum', 'pvivax', 'pknowlesi', 'pberghei', 'tgondii']
COL = "query,target,qstart,qend,tstart,tend,qlen,tlen,pident,evalue,bits,qaln,taln".split(",")
STRUCT_MAP = {'internal_insertion': 'insertion', 'internal_variable': 'opposite_segment',
              'terminal_extension': 'terminal', 'terminal_variable': 'terminal'}


def rd(p):
    with open(p, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def flip(h):
    """Swap query/target of an m8 hit dict."""
    return {'query': h['target'], 'target': h['query'], 'qstart': h['tstart'], 'tstart': h['qstart'],
            'qaln': h['taln'], 'taln': h['qaln']}


def call(h, s, e):
    amap = aligned_map(h)
    if not amap:
        return 'no_alignment', None, None
    L, R = anchors(amap, s, e)
    if L is None or R is None:
        return 'terminal', None, None
    g = amap[R] - amap[L] - 1
    a_span = R - L - 1
    return ('insertion' if g <= max(5, 0.25 * a_span) else 'opposite_segment'), g, a_span


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rbh', required=True)
    ap.add_argument('--lcr', required=True)
    ap.add_argument('--placement', required=True)
    ap.add_argument('--status', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    has_model = {(r['species'], r['gene']) for r in rd(a.status) if r['file']}
    regions = collections.defaultdict(list)
    for sp in PILOT:
        for r in rd(os.path.join(a.lcr, f'{sp}_lcr.tsv')):
            regions[(sp, r['gene'])].append((int(r['start']) - 1, int(r['end']), r['lcr_class']))
    struct = {}
    for r in rd(os.path.join(a.placement, 'region_partner.tsv')):
        if r['location'] in STRUCT_MAP:
            struct[(r['a_species'], r['a_gene'], r['b_species'], int(r['region_start']), int(r['region_end']))] = STRUCT_MAP[r['location']]
    struct_idx = collections.defaultdict(list)
    for (asp, ag, bsp, s, e), v in struct.items():
        struct_idx[(asp, ag, bsp)].append((s, e, v))

    rows = []
    for i, A in enumerate(PILOT):
        for B in PILOT[i + 1:]:
            p = os.path.join(a.rbh, f'rbh_{A}__{B}.tsv')
            if not os.path.exists(p):
                continue
            for line in open(p):
                f = line.rstrip('\n').split('\t')
                if len(f) < 13:
                    continue
                h = dict(zip(COL, f))
                for (sa, ga, sb, gb, hh) in ((A, h['query'], B, h['target'], h), (B, h['target'], A, h['query'], flip(h))):
                    for s, e, cls in regions.get((sa, ga), []):
                        c, g, span = call(hh, s, e)
                        sc = None
                        for s2, e2, v in struct_idx.get((sa, ga, sb), []):
                            if abs(s2 - (s + 1)) <= 3 and abs(e2 - e) <= 3:
                                sc = v
                                break
                        rows.append(dict(a_species=sa, a_gene=ga, b_species=sb, b_gene=gb, start=s + 1, end=e,
                                         lcr_class=cls, seq_call=c, b_seg=g, a_span=span, struct_call=sc or '',
                                         both_models=(sa, ga) in has_model and (sb, gb) in has_model))
    # validation
    conf = collections.Counter((r['struct_call'], r['seq_call']) for r in rows if r['struct_call'])
    val = [r for r in rows if r['struct_call'] and r['seq_call'] in ('insertion', 'opposite_segment', 'terminal')]
    agree = sum(r['struct_call'] == r['seq_call'] for r in val)
    internal = [r for r in val if r['struct_call'] in ('insertion', 'opposite_segment') and r['seq_call'] in ('insertion', 'opposite_segment')]
    ins_s = [r for r in internal if r['struct_call'] == 'insertion']
    ins_q = [r for r in internal if r['seq_call'] == 'insertion']
    tp = sum(1 for r in internal if r['struct_call'] == 'insertion' and r['seq_call'] == 'insertion')

    def dist(sub):
        c = collections.Counter(r['seq_call'] for r in sub)
        n = sum(c[k] for k in ('insertion', 'opposite_segment', 'terminal'))
        return dict(n=n, **{k: round(c[k] / n, 3) for k in ('insertion', 'opposite_segment', 'terminal')}) if n else {}
    nomodel = [r for r in rows if not r['both_models']]
    withmodel = [r for r in rows if r['both_models']]
    by_class = {}
    for cls in sorted({r['lcr_class'] for r in rows}):
        by_class[cls] = dict(with_models=dist([r for r in withmodel if r['lcr_class'] == cls]),
                             without_models=dist([r for r in nomodel if r['lcr_class'] == cls]))
    res = dict(
        n_rows=len(rows),
        validation=dict(n=len(val), overall_agreement=round(agree / len(val), 3),
                        internal_n=len(internal),
                        insertion_recall=round(tp / len(ins_s), 3) if ins_s else None,
                        insertion_precision=round(tp / len(ins_q), 3) if ins_q else None,
                        confusion={f'{k[0]}->{k[1]}': v for k, v in sorted(conf.items())}),
        seq_placement_with_models=dist(withmodel),
        seq_placement_without_models=dist(nomodel),
        by_class=by_class)
    json.dump(res, open(a.out + '.json', 'w'), indent=1)
    with open(a.out + '.tsv', 'w') as fh:
        ks = list(rows[0].keys())
        fh.write('\t'.join(ks) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[k]) for k in ks) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
