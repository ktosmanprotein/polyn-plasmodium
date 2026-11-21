"""Structure helpers for the poly-N study: read AlphaFold PDB(.gz) models,
map coordinates between annotation and model sequences, and a light
CA-geometry secondary-structure assignment. Standard library only.
"""
from __future__ import annotations

import difflib
import gzip
import math

THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
         "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
         "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "MSE": "M"}


def read_model(path):
    """Return dict(seq, ca=[(x,y,z)], plddt=[float]) from an AlphaFold PDB(.gz).

    One entry per residue (first chain, CA atoms), in model order. AlphaFold
    stores per-residue pLDDT in the B-factor column."""
    op = gzip.open if str(path).endswith(".gz") else open
    seq, ca, pl = [], [], []
    chain = None
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            ch = line[21]
            if chain is None:
                chain = ch
            elif ch != chain:
                break
            seq.append(THREE.get(line[17:20].strip(), "X"))
            ca.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            pl.append(float(line[60:66]))
    return {"seq": "".join(seq), "ca": ca, "plddt": pl}


def coord_map(src, dst):
    """Map 0-based positions of sequence `src` onto `dst`.

    Returns (mapping list of length len(src) with dst index or None, identity
    over src). Identical sequences map 1:1; otherwise difflib matching blocks
    (sufficient for the small annotation differences between VEuPathDB and
    UniProt; heavily divergent pairs show up as low identity and are dropped
    downstream)."""
    if src == dst:
        return list(range(len(src))), 1.0
    m = [None] * len(src)
    sm = difflib.SequenceMatcher(None, src, dst, autojunk=False)
    same = 0
    for a, b, n in sm.get_matching_blocks():
        for k in range(n):
            m[a + k] = b + k
        same += n
    return m, same / max(len(src), 1)


def dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def ca_sse(ca):
    """Secondary structure from CA geometry (P-SEA-like): 'H' helix, 'E'
    strand, 'C' other. Uses i..i+2 / i..i+3 / i..i+4 CA distances.

    Helix: d(i,i+2)~5.4, d(i,i+3)~5.0, d(i,i+4)~6.2. Strand: d(i,i+2)~6.7,
    d(i,i+3)~9.9, d(i,i+4)~12.4. A residue is assigned when a window of 5
    residues starting at it meets all three, then helix runs < 4 and strand
    runs < 3 are cleared. Approximate by design; used only to say whether a
    position sits inside a regular element or a loop."""
    n = len(ca)
    s = ["C"] * n
    for i in range(n - 4):
        d2, d3, d4 = dist(ca[i], ca[i + 2]), dist(ca[i], ca[i + 3]), dist(ca[i], ca[i + 4])
        if abs(d2 - 5.45) < 0.6 and abs(d3 - 5.1) < 0.75 and abs(d4 - 6.2) < 0.85:
            for k in range(i, i + 5):
                s[k] = "H"
        elif abs(d2 - 6.4) < 0.8 and abs(d3 - 9.9) < 1.3 and abs(d4 - 12.4) < 1.6:
            for k in range(i, i + 5):
                if s[k] != "H":
                    s[k] = "E"
    # a strand residue must pair with a residue elsewhere in a sheet:
    # a non-local CA (|i-j| > 3) within 5.5 A. Extended coil (common in
    # disordered, low-pLDDT regions) has no partner and is reset to coil.
    grid = {}
    for i, p in enumerate(ca):
        grid.setdefault((int(p[0] // 6), int(p[1] // 6), int(p[2] // 6)), []).append(i)
    for i in range(n):
        if s[i] != "E":
            continue
        gx, gy, gz = (int(c // 6) for c in ca[i])
        paired = False
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get((gx + dx, gy + dy, gz + dz), ()):
                        if abs(i - j) > 3 and dist(ca[i], ca[j]) < 5.5:
                            paired = True
                            break
                    if paired:
                        break
                if paired:
                    break
            if paired:
                break
        if not paired:
            s[i] = "C"
    out = "".join(s)
    # clear runs that are too short to be real elements
    res = list(out)
    i = 0
    while i < n:
        j = i
        while j < n and out[j] == out[i]:
            j += 1
        if (out[i] == "H" and j - i < 4) or (out[i] == "E" and j - i < 3):
            for k in range(i, j):
                res[k] = "C"
        i = j
    return "".join(res)


def aligned_pairs(qaln, taln, qstart, tstart):
    """From Foldseek qaln/taln strings (1-based starts), yield (qpos, tpos)
    0-based for every column where both residues are present."""
    qi, ti = qstart - 1, tstart - 1
    for a, b in zip(qaln, taln):
        if a != "-" and b != "-":
            yield qi, ti
        if a != "-":
            qi += 1
        if b != "-":
            ti += 1
