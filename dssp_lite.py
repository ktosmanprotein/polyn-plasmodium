"""Minimal DSSP (Kabsch & Sander 1983) in numpy, 3-state output.

Backbone N, CA, C, O from an AlphaFold PDB(.gz); amide H placed DSSP-style
(N + unit(C_prev - O_prev)). H-bond if electrostatic energy
E = 0.084 * 332 * (1/rON + 1/rCH - 1/rOH - 1/rCN) < -0.5 kcal/mol (O,C of residue i
as acceptor; N,H of residue j as donor).
  H (alpha helix): two consecutive 4-turns at i-1 and i -> residues i..i+3 = H
  G (3-10) and I (pi) handled the same way with 3- and 5-turns, only where not H
  E/B: parallel / antiparallel bridges by the DSSP patterns
3-state: H,G,I -> 'H'; E,B -> 'E'; rest 'C'.
"""
import gzip
import numpy as np


def read_backbone(path):
    op = gzip.open if str(path).endswith('.gz') else open
    res = {}
    order = []
    chain = None
    pl = {}
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith('ENDMDL'):
                break
            if not line.startswith('ATOM'):
                continue
            ch = line[21]
            if chain is None:
                chain = ch
            elif ch != chain:
                break
            name = line[12:16].strip()
            if name not in ('N', 'CA', 'C', 'O'):
                continue
            key = (int(line[22:26]), line[26])
            if key not in res:
                res[key] = {}
                order.append(key)
            res[key][name] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if name == 'CA':
                pl[key] = float(line[60:66])
    keep = [k for k in order if all(a in res[k] for a in ('N', 'CA', 'C', 'O'))]
    X = np.array([[res[k][a] for a in ('N', 'CA', 'C', 'O')] for k in keep], dtype=np.float64)
    return X, np.array([pl[k] for k in keep])


def hbond_map(X):
    N, CA, C, O = X[:, 0], X[:, 1], X[:, 2], X[:, 3]
    n = len(X)
    H = N.copy()
    v = C[:-1] - O[:-1]
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    H[1:] = N[1:] + v
    # acceptor i (C,O) -> donor j (N,H); computed in row blocks to bound memory
    def d(A, B):
        return np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)) + 1e-6
    hb = np.zeros((n, n), bool)
    for s in range(0, n, 400):
        e = min(n, s + 400)
        Ob, Cb = O[s:e], C[s:e]
        # only pairs with O..N < 5.2 A can pass the cutoff in practice; full formula kept
        E = 0.084 * 332 * (1 / d(Ob, N) + 1 / d(Cb, H) - 1 / d(Ob, H) - 1 / d(Cb, N))
        hb[s:e] = E < -0.5
    idx = np.arange(n)
    hb[np.abs(idx[:, None] - idx[None, :]) < 2] = False
    hb[:, 0] = False          # first residue has no H
    return hb                 # hb[i, j]: O(i) accepts from N-H(j)


def assign(X):
    n = len(X)
    if n < 5:
        return 'C' * n
    hb = hbond_map(X)
    ss = np.array(['C'] * n)

    def turns(k):
        t = np.zeros(n, bool)
        i = np.arange(n - k)
        t[i] = hb[i, i + k]
        return t

    for k, lab in ((4, 'H'), (3, 'G'), (5, 'I')):
        t = turns(k)
        for i in range(1, n - k):
            if t[i - 1] and t[i]:
                for m in range(i, i + k):
                    if lab == 'H' or ss[m] == 'C':
                        ss[m] = lab
    # bridges
    T = hb.T
    P = np.zeros((n + 2, n + 2), bool)
    Q = np.zeros((n + 2, n + 2), bool)
    P[1:-1, 1:-1] = hb            # P[i+1, j+1] = hb[i, j]
    Q[1:-1, 1:-1] = T
    I = slice(1, n + 1)           # padded index of residue i
    Im, Ip = slice(0, n), slice(2, n + 2)
    par1 = P[Im, I] & Q[Ip, I]            # hb[i-1, j] & hb[j, i+1]
    par2 = Q[I, Im] & P[I, Ip]            # hb[j-1, i] & hb[i, j+1]
    anti1 = hb & T                        # hb[i, j] & hb[j, i]
    anti2 = P[Im, Ip] & Q[Ip, Im]         # hb[i-1, j+1] & hb[j-1, i+1]
    B = par1 | par2 | anti1 | anti2
    idx = np.arange(n)
    B[np.abs(idx[:, None] - idx[None, :]) < 3] = False
    B[0, :] = B[-1, :] = False
    bridge = B.any(1)
    for i in range(n):
        if bridge[i] and ss[i] == 'C':
            ss[i] = 'E'
    out = ''.join('H' if c in 'HGI' else ('E' if c == 'E' else 'C') for c in ss)
    return out
