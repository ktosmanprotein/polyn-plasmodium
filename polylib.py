"""Shared logic for the poly-N / low-complexity insertion study.

Standard library only, so it runs under any python3 >= 3.8.
"""
from __future__ import annotations

import csv
import gzip
import math
import os
import re
from collections import Counter

# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_fasta(path):
    """Yield (header_without_>, sequence) pairs."""
    header, chunks = None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".part"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# VEuPathDB headers
# ---------------------------------------------------------------------------

def parse_veupath_header(header):
    """Parse '>ID | key=value | key=value ...' into a dict.

    VEuPathDB protein headers look like
      PF3D7_0100100.1-p1 | transcript=PF3D7_0100100.1 | gene=PF3D7_0100100 | ...
    The first token is kept as 'id'. If no gene= field exists, the gene is
    taken from the id by stripping a transcript suffix (.1, .1-p1).
    """
    parts = [p.strip() for p in header.split("|")]
    out = {"id": parts[0].split()[0] if parts and parts[0] else ""}
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            out[k.strip()] = v.strip()
    if "gene" not in out:
        out["gene"] = re.sub(r"(\.\d+)?(-p\d+)?$", "", out["id"])
    if "transcript" not in out:
        out["transcript"] = re.sub(r"-p\d+$", "", out["id"])
    return out


def is_pseudo(meta):
    v = str(meta.get("is_pseudo", "")).lower()
    return v in ("true", "1", "yes")


def primary_isoforms(records):
    """Keep one record per gene: the lexically first transcript (.1 before .2).

    records: iterable of (meta, seq). Returns dict gene -> (meta, seq).
    """
    best = {}
    for meta, seq in records:
        g = meta["gene"]
        if g not in best or meta["transcript"] < best[g][0]["transcript"]:
            best[g] = (meta, seq)
    return best


# ---------------------------------------------------------------------------
# Homopolymer runs
# ---------------------------------------------------------------------------

def homopolymer_runs(seq, residue, min_len):
    """Maximal runs of `residue` of length >= min_len, as (start, end) with
    0-based start and exclusive end."""
    pat = re.compile(re.escape(residue) + "{%d,}" % min_len)
    return [(m.start(), m.end()) for m in pat.finditer(seq)]


def expected_maximal_runs(length, p, k):
    """Expected number of maximal runs of length >= k in an i.i.d. sequence
    of `length` positions where each position is the residue with prob p.

    A run of >= k starts at position 0 with prob p^k, and at any later
    position i <= length-k with prob (1-p) p^k (previous position differs).
    """
    if length < k or p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    return (p ** k) * (1.0 + (length - k) * (1.0 - p))


# ---------------------------------------------------------------------------
# Low-complexity regions (SEG-like)
# ---------------------------------------------------------------------------

def _entropy_from_counts(counts, n):
    h = 0.0
    for c in counts.values():
        if c:
            q = c / n
            h -= q * math.log2(q)
    return h


def window_entropies(seq, window):
    """Shannon entropy (bits) of every length-`window` window, by start."""
    n = len(seq)
    if n < window:
        return []
    counts = Counter(seq[:window])
    out = [_entropy_from_counts(counts, window)]
    for i in range(1, n - window + 1):
        old, new = seq[i - 1], seq[i + window - 1]
        if old != new:
            counts[old] -= 1
            if counts[old] == 0:
                del counts[old]
            counts[new] += 1
        out.append(_entropy_from_counts(counts, window))
    return out


def lcr_regions(seq, window=12, trigger=2.2, extend=2.5, min_len=0):
    """SEG-like low-complexity regions.

    Windows with entropy <= `extend` are 'low'. A contiguous block of low
    windows becomes a region if it contains at least one window at or below
    `trigger`; the region spans from the first window's start to the last
    window's end. This follows SEG's two-threshold logic (trigger, then
    extend) without SEG's final boundary optimisation, and is stated as such
    in Methods. Returns sorted, non-overlapping (start, end) pairs.
    """
    ent = window_entropies(seq, window)
    regions = []
    i, m = 0, len(ent)
    while i < m:
        if ent[i] <= extend:
            j = i
            triggered = False
            while j < m and ent[j] <= extend:
                if ent[j] <= trigger:
                    triggered = True
                j += 1
            if triggered:
                regions.append((i, (j - 1) + window))
            i = j
        else:
            i += 1
    # windows overlap, so blocks separated by one high window can still
    # overlap in residue space: merge.
    merged = []
    for s, e in regions:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return [(s, e) for s, e in merged if e - s >= min_len]


def classify_segment(segment, rules):
    """Composition class of a segment. `rules` is an ordered list of
    (label, residues, min_fraction); the first rule met wins."""
    n = len(segment)
    if n == 0:
        return "empty"
    c = Counter(segment)
    for label, residues, frac in rules:
        if sum(c[r] for r in residues) / n >= frac:
            return label
    return "mixed"


def top_composition(segment, k=3):
    n = len(segment)
    c = Counter(segment).most_common(k)
    return ",".join(f"{aa}:{cnt / n:.2f}" for aa, cnt in c)


DEFAULT_CLASS_RULES = [
    ("N-rich", "N", 0.40),
    ("Q-rich", "Q", 0.40),
    ("GS-rich", "GS", 0.50),
    ("charged-KED", "KED", 0.50),
]


# ---------------------------------------------------------------------------
# Mutation-bias null from coding sequence
# ---------------------------------------------------------------------------

#: Codon prefixes whose third position is fourfold degenerate.
FOURFOLD_PREFIXES = ("GC", "CG", "GG", "CT", "CC", "AC", "TC", "GT")
STOPS = ("TAA", "TAG", "TGA")


def fourfold_composition(cds_iter):
    """Nucleotide composition at fourfold-degenerate third positions across
    an iterable of CDS strings. Returns dict A/C/G/T -> frequency, and the
    number of sites counted."""
    counts = Counter()
    for cds in cds_iter:
        cds = cds.upper()
        if len(cds) % 3:
            continue
        for i in range(0, len(cds) - 2, 3):
            codon = cds[i:i + 3]
            if codon[:2] in FOURFOLD_PREFIXES and codon[2] in "ACGT":
                counts[codon[2]] += 1
    total = sum(counts[b] for b in "ACGT")
    if total == 0:
        return {b: 0.25 for b in "ACGT"}, 0
    return {b: counts[b] / total for b in "ACGT"}, total


def residue_prob_under_composition(freqs, codons):
    """Probability that a random sense codon drawn from i.i.d. nucleotides
    with frequencies `freqs` encodes one of `codons`."""
    def pc(c):
        return freqs[c[0]] * freqs[c[1]] * freqs[c[2]]
    p_stop = sum(pc(s) for s in STOPS)
    return sum(pc(c) for c in codons) / (1.0 - p_stop)


CODONS = {"N": ("AAT", "AAC"), "K": ("AAA", "AAG"), "Q": ("CAA", "CAG")}
