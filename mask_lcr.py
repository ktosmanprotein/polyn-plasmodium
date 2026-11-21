#!/usr/bin/env python3
"""Step 2b (prep) — write LCR-masked primary-isoform proteomes for the
homology search.

Each fLPS region from step 2a is replaced by X. Masking keeps homology calls
from being driven by shared composition (two unrelated N-rich proteins look
similar to an aligner) and keeps the search tractable: unmasked, a single
P. falciparum vs P. berghei search ran >20 min without finishing on 2 cores.
Orthology should rest on the folded parts of the proteins in any case, since
the inserts are what the study measures.

Output headers are gene ids, one record per gene (primary isoform).
"""
import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polylib as P  # noqa: E402


def mask(seq, regions):
    s = list(seq)
    for a, b in regions:
        for i in range(max(a, 0), min(b, len(s))):
            s[i] = "X"
    return "".join(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--seqs", required=True)
    ap.add_argument("--lcr-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-unmasked", type=int, default=30,
                    help="drop proteins with fewer unmasked residues than this")
    a = ap.parse_args()
    from call_tracts import load_proteins  # noqa: E402

    os.makedirs(a.out, exist_ok=True)
    for sp in P.read_tsv(a.species):
        key = sp["key"]
        prot = os.path.join(a.seqs, f"{key}_proteins.fasta")
        lcr = os.path.join(a.lcr_dir, f"{key}_lcr.tsv")
        if not (os.path.exists(prot) and os.path.exists(lcr)):
            print(f"[skip] {key}")
            continue
        proteins, _ = load_proteins(prot)
        regs = defaultdict(list)
        for r in P.read_tsv(lcr):
            regs[r["gene"]].append((int(r["start"]) - 1, int(r["end"])))
        out = os.path.join(a.out, f"{key}_masked.fasta")
        kept = dropped = masked_res = total = 0
        with open(out + ".part", "w") as fh:
            for g, (_, seq) in sorted(proteins.items()):
                m = mask(seq, regs.get(g, []))
                total += len(seq)
                masked_res += m.count("X") - seq.count("X")
                if len(m) - m.count("X") < a.min_unmasked:
                    dropped += 1
                    continue
                fh.write(f">{g}\n")
                for i in range(0, len(m), 80):
                    fh.write(m[i:i + 80] + "\n")
                kept += 1
        os.replace(out + ".part", out)
        print(f"[done] {key}: {kept:,} proteins written, {dropped} dropped "
              f"(< {a.min_unmasked} unmasked residues); {100 * masked_res / total:.1f}% residues masked")
    return 0


if __name__ == "__main__":
    sys.exit(main())
