#!/usr/bin/env python3
"""Item 1 — gain/loss reconstruction of Pf poly-N tracts across 13 Plasmodium species.

Input: step-6 dating table (pf_tract_dating[_lenient].tsv). Per tract and species:
  present (1), present_weak (1 in lenient mode, 0 in strict), absent / absent_noseg (0),
  no_ortholog (missing). Outside Laverania 'present' = N-rich at the homologous segment.

Three rootings of the contested Laverania placement (Laverania internal tree after Otto 2018):
  T1 mammal  : ((Lav, (primate, rodent)), avian)       Laverania inside mammalian Plasmodium
  T2 avian   : ((Lav, avian), (primate, rodent))       Laverania sister to avian parasites
  T3 basal   : (Lav, ((primate, rodent), avian))       Laverania sister to all others

For each tree:
  * Fitch parsimony: steps per tract, total steps (which topology the tract patterns favour).
  * ML two-state Mk model (gain g, loss l; one shared Q, free branch lengths), likelihood
    conditioned on presence in Pf (all tracts are ascertained in Pf). Marginal posterior of
    presence at each internal node; expected gains/losses per branch by posterior smoothing.
  Readouts: P(present) at the Laverania crown (Lav) and its parent (stem, LavP); for tracts
  that are N-rich at the same site outside Laverania, retention (present at LavP) vs
  independent gains.
"""
import argparse, json, math, collections
import numpy as np
from scipy.optimize import minimize

SP = ['Pf', 'Pprf', 'Pr', 'Pbilc', 'Pblac', 'Padl', 'Pgab', 'Pv', 'Pk', 'Pb', 'Py', 'Pgal', 'Prel']
LAV = "(((((Pf,Pprf)n_FP,Pr)n_FPR,Pbilc)n_FPRB,Pblac)CladeB,(Padl,Pgab)CladeA)Lav"
PRI, ROD, AVI = "(Pv,Pk)Primate", "(Pb,Py)Rodent", "(Pgal,Prel)Avian"
TREES = {
    'T1_mammal': f"(({LAV},({PRI},{ROD})Mammal_nonLav)LavP,{AVI})Root",
    'T2_avian': f"(({LAV},{AVI})LavP,({PRI},{ROD})Mammal_nonLav)Root",
    'T3_basal': f"({LAV},(({PRI},{ROD})Mammal_nonLav,{AVI})Others)LavP",
}


def parse(nwk):
    """Return nodes list (post-order) of dicts {name, children}."""
    pos = 0
    nodes = []

    def node():
        nonlocal pos
        ch = []
        if nwk[pos] == '(':
            pos += 1
            while True:
                ch.append(node())
                if nwk[pos] == ',':
                    pos += 1
                    continue
                if nwk[pos] == ')':
                    pos += 1
                    break
        start = pos
        while pos < len(nwk) and nwk[pos] not in ',();':
            pos += 1
        name = nwk[start:pos]
        idx = len(nodes)
        nodes.append({'name': name, 'children': ch})
        return idx

    root = node()
    return nodes, root


def load(path, lenient):
    rows = []
    with open(path) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for line in f:
            d = dict(zip(hdr, line.rstrip('\n').split('\t')))
            st = {'Pf': 1}
            for s in SP[1:]:
                v = d.get(s, '')
                if v == 'present':
                    st[s] = 1
                elif v == 'present_weak':
                    st[s] = 1 if lenient else 0
                elif v in ('absent', 'absent_noseg'):
                    st[s] = 0
                else:
                    st[s] = None
            rows.append({'gene': d['gene'], 'start': d['start'], 'origin': d['origin'],
                         'len': int(d['length']), 'st': st})
    return rows


def fitch(nodes, root, st):
    sets = {}
    steps = 0
    for i, n in enumerate(nodes):
        if not n['children']:
            v = st.get(n['name'])
            sets[i] = {0, 1} if v is None else {v}
        else:
            a, b = (sets[c] for c in n['children'])
            if a & b:
                sets[i] = a & b
            else:
                sets[i] = a | b
                steps += 1
    return steps, sets


def leaf_partials(nodes, rows):
    P = {}
    for i, n in enumerate(nodes):
        if not n['children']:
            m = np.ones((len(rows), 2))
            for k, r in enumerate(rows):
                v = r['st'].get(n['name'])
                if v == 0:
                    m[k] = [1, 0]
                elif v == 1:
                    m[k] = [0, 1]
            P[i] = m
    return P


def pmat(g, l, t):
    s = g + l
    e = math.exp(-s * t)
    pi1 = g / s
    pi0 = l / s
    return np.array([[pi0 + pi1 * e, pi1 - pi1 * e],
                     [pi0 - pi0 * e, pi1 + pi0 * e]])


class Model:
    def __init__(self, nodes, root, rows):
        self.nodes, self.root, self.rows = nodes, root, rows
        self.leafP = leaf_partials(nodes, rows)
        self.branches = [i for i in range(len(nodes)) if i != root]
        self.parent = {}
        for i, n in enumerate(nodes):
            for c in n['children']:
                self.parent[c] = i
        self.pf = next(i for i, n in enumerate(nodes) if n['name'] == 'Pf')
        # path from root to Pf for the ascertainment term
        p, path = self.pf, []
        while p != root:
            path.append(p)
            p = self.parent[p]
        self.pf_path = path[::-1]

    def unpack(self, x):
        g, l = math.exp(x[0]), math.exp(x[1])
        s = g + l
        g, l = g / s, l / s                   # rates normalised; time absorbed into t
        t = {b: math.exp(x[2 + k]) for k, b in enumerate(self.branches)}
        return g, l, t

    def up(self, g, l, t):
        L = {}
        for i, n in enumerate(self.nodes):
            if not n['children']:
                L[i] = self.leafP[i]
            else:
                acc = np.ones((len(self.rows), 2))
                for c in n['children']:
                    acc *= L[c] @ pmat(g, l, t[c]).T
                L[i] = acc
        return L

    def nll(self, x):
        g, l, t = self.unpack(x)
        pi = np.array([l, g])
        L = self.up(g, l, t)
        lik = L[self.root] @ pi
        # P(Pf = 1)
        v = pi.copy()
        for b in self.pf_path:
            v = v @ pmat(g, l, t[b])
        asc = v[1]
        return -(np.log(np.maximum(lik, 1e-300)).sum() - len(self.rows) * math.log(asc))

    def fit(self):
        x0 = np.array([math.log(0.3), math.log(0.7)] + [math.log(0.2)] * len(self.branches))
        res = minimize(self.nll, x0, method='L-BFGS-B',
                       bounds=[(-8, 3), (-8, 3)] + [(-9, 3)] * len(self.branches))
        self.x = res.x
        self.loglik = -res.fun
        return res

    def marginals(self):
        g, l, t = self.unpack(self.x)
        pi = np.array([l, g])
        L = self.up(g, l, t)
        n = len(self.rows)
        # downward (outside) messages
        O = {self.root: np.tile(pi, (n, 1))}
        order = list(range(len(self.nodes)))[::-1]   # pre-order (reverse post-order)
        for i in order:
            ch = self.nodes[i]['children']
            for c in ch:
                sib = np.ones((n, 2))
                for c2 in ch:
                    if c2 != c:
                        sib *= L[c2] @ pmat(g, l, t[c2]).T
                O[c] = (O[i] * sib) @ pmat(g, l, t[c])
        post = {}
        for i in range(len(self.nodes)):
            m = O[i] * L[i]
            post[i] = m[:, 1] / m.sum(1)
        # expected gains / losses per branch (joint parent-child posterior)
        exp_gain, exp_loss = {}, {}
        for c in self.branches:
            p = self.parent[c]
            sib = np.ones((n, 2))
            for c2 in self.nodes[p]['children']:
                if c2 != c:
                    sib *= L[c2] @ pmat(g, l, t[c2]).T
            outside = O[p] * sib                       # n x 2 over parent state
            P = pmat(g, l, t[c])
            joint = outside[:, :, None] * P[None] * L[c][:, None, :]
            joint /= joint.sum((1, 2), keepdims=True)
            exp_gain[c] = joint[:, 0, 1]
            exp_loss[c] = joint[:, 1, 0]
        return g, l, t, post, exp_gain, exp_loss


def outside_present(r):
    return any(r['st'].get(s) == 1 for s in ['Pv', 'Pk', 'Pb', 'Py', 'Pgal', 'Prel'])


def run(rows, label):
    out = {}
    for tn, nwk in TREES.items():
        nodes, root = parse(nwk)
        name2i = {n['name']: i for i, n in enumerate(nodes)}
        steps = []
        fitch_lavp = collections.Counter()
        for r in rows:
            s, sets = fitch(nodes, root, r['st'])
            steps.append(s)
            fitch_lavp[str(sorted(sets[name2i['LavP']]))] += 1
        m = Model(nodes, root, rows)
        m.fit()
        g, l, t, post, eg, el = m.marginals()
        lav, lavp = post[name2i['Lav']], post[name2i['LavP']]
        allLav = np.array([r['origin'].startswith('all Laverania') for r in rows])
        shared = np.array([outside_present(r) for r in rows])
        nonlav_branches = [b for b in m.branches if nodes[b]['name'] in
                           ('Pv', 'Pk', 'Pb', 'Py', 'Pgal', 'Prel', 'Primate', 'Rodent', 'Avian',
                            'Mammal_nonLav', 'Others')]
        gains_out = sum(eg[b] for b in nonlav_branches)
        res = {
            'loglik': round(m.loglik, 2), 'gain_rate': round(g, 3), 'loss_rate': round(l, 3),
            'branch_lengths': {nodes[b]['name']: round(t[b], 3) for b in m.branches},
            'fitch_total_steps': int(sum(steps)),
            'fitch_LavP_state_counts': dict(fitch_lavp),
            'n_tracts': len(rows),
            'n_allLav': int(allLav.sum()), 'n_allLav_shared_outside': int((allLav & shared).sum()),
            'mean_P_present_Lav_crown_allLav': round(float(lav[allLav].mean()), 3),
            'n_allLav_present_at_LavP_gt0.5': int((lavp[allLav] > 0.5).sum()),
            'shared_outside_allLav': {
                'n': int((allLav & shared).sum()),
                'retained_from_LavP_gt0.5': int((lavp[allLav & shared] > 0.5).sum()),
                'mean_P_LavP': round(float(lavp[allLav & shared].mean()), 3),
                'expected_independent_gains_outside_Lav': round(float(gains_out[allLav & shared].sum()), 1),
            },
            'not_shared_allLav_mean_P_LavP': round(float(lavp[allLav & ~shared].mean()), 3),
            'expected_gains_by_branch': {nodes[b]['name']: round(float(eg[b].sum()), 1) for b in m.branches},
            'expected_losses_by_branch': {nodes[b]['name']: round(float(el[b].sum()), 1) for b in m.branches},
            'node_mean_P_present': {nodes[i]['name']: round(float(post[i].mean()), 3)
                                    for i in range(len(nodes)) if nodes[i]['children']},
        }
        out[tn] = res
        out.setdefault('_per_tract', {})[tn] = [round(float(x), 3) for x in lavp]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--strict', required=True)
    ap.add_argument('--lenient', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    result = {}
    for label, path, len_ in [('strict', a.strict, False), ('lenient', a.lenient, True)]:
        rows = load(path, len_)
        rows = [r for r in rows if r['origin'] != 'uninformative']
        r = run(rows, label)
        per = r.pop('_per_tract')
        result[label] = r
        with open(f'{a.out}_{label}_per_tract.tsv', 'w') as f:
            f.write('gene\tstart\torigin\tlen\t' + '\t'.join(f'P_LavP_{k}' for k in per) + '\n')
            for i, row in enumerate(rows):
                f.write(f"{row['gene']}\t{row['start']}\t{row['origin']}\t{row['len']}\t" +
                        '\t'.join(str(per[k][i]) for k in per) + '\n')
    json.dump(result, open(a.out + '.json', 'w'), indent=1)
    for label in result:
        for tn, r in result[label].items():
            print(label, tn, 'logL', r['loglik'], 'fitch', r['fitch_total_steps'],
                  'g/l', r['gain_rate'], r['loss_rate'], 'shared', r['shared_outside_allLav'])


if __name__ == '__main__':
    main()
