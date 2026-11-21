#!/usr/bin/env python3
"""Step 2b — one-to-one orthogroups from all-pairs reciprocal best hits.

Every pair of species is searched with MMseqs2 easy-rbh (E <= 1e-5, low-
complexity masking on, which is MMseqs2's default). Reciprocal best hits are
edges in a graph of genes; connected components are candidate orthogroups.

A component is kept as a 1:1 orthogroup only if it has at most one gene per
species. Components that pull in two genes from one species (paralogues
chained through RBH edges) are reported as 'conflict' and left out, rather
than resolved by a rule that would be hard to defend.

Output
  orthogroups.tsv      og_id, status, n_species, plasmodium_complete, and one
                       gene column plus one UniProt column per species
  orthogroup_summary   counts by status and species coverage
"""
import argparse
import glob
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polylib as P  # noqa: E402


class DSU:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


def id_to_gene(fasta):
    out = {}
    for h, _ in P.read_fasta(fasta):
        m = P.parse_veupath_header(h)
        out[m["id"]] = m["gene"]
    return out


def build(species, rbh_files, id2gene):
    """species: ordered keys. rbh_files: {(a, b): path}. id2gene: {key: {id: gene}}.
    Returns list of orthogroup dicts."""
    dsu = DSU()
    edges = 0
    for (a, b), path in rbh_files.items():
        with open(path) as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) < 2:
                    continue
                ga = id2gene[a].get(p[0])
                gb = id2gene[b].get(p[1])
                if ga and gb:
                    dsu.union((a, ga), (b, gb))
                    edges += 1
    comps = defaultdict(list)
    for node in list(dsu.p):
        comps[dsu.find(node)].append(node)

    plas = [s for s in species if s.startswith("p")]
    groups = []
    for nodes in comps.values():
        by_sp = defaultdict(list)
        for sp, g in nodes:
            by_sp[sp].append(g)
        status = "1to1" if all(len(v) == 1 for v in by_sp.values()) else "conflict"
        groups.append({
            "status": status,
            "members": {sp: sorted(v) for sp, v in by_sp.items()},
            "n_species": len(by_sp),
            "plasmodium_complete": int(all(s in by_sp for s in plas)),
        })
    # stable ids: sort by the P. falciparum gene (or first available gene)
    def sort_key(gp):
        for s in species:
            if s in gp["members"]:
                return (species.index(s), gp["members"][s][0])
        return (99, "")
    groups.sort(key=sort_key)
    for i, gp in enumerate(groups, 1):
        gp["og_id"] = f"OG{i:05d}"
    return groups, edges


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--seqs", required=True, help="folder with the searched FASTA files")
    ap.add_argument("--fasta-suffix", default="_masked.fasta",
                    help="file name after the species key (default: step 2b masked files)")
    ap.add_argument("--rbh-dir", required=True, help="folder with rbh_<a>__<b>.m8")
    ap.add_argument("--apico-proteomes", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    from call_tracts import uniprot_map  # noqa: E402

    species = [r["key"] for r in P.read_tsv(a.species)]
    id2gene = {s: id_to_gene(os.path.join(a.seqs, f"{s}{a.fasta_suffix}")) for s in species}
    rbh = {}
    for path in glob.glob(os.path.join(a.rbh_dir, "rbh_*__*.m8")):
        x, y = os.path.basename(path)[4:-3].split("__")
        if x in species and y in species:
            rbh[(x, y)] = path
    expected = len(species) * (len(species) - 1) // 2
    print(f"{len(rbh)} of {expected} species pairs present")
    if len(rbh) < expected:
        print("WARNING: missing pairs make orthogroups incomplete")

    groups, edges = build(species, rbh, id2gene)
    uni = {s: uniprot_map(a.apico_proteomes, s) for s in species}

    fields = ["og_id", "status", "n_species", "plasmodium_complete"]
    fields += [f"{s}_gene" for s in species] + [f"{s}_uniprot" for s in species]
    rows = []
    for gp in groups:
        r = {k: gp[k] for k in ("og_id", "status", "n_species", "plasmodium_complete")}
        for s in species:
            genes = gp["members"].get(s, [])
            r[f"{s}_gene"] = ",".join(genes)
            r[f"{s}_uniprot"] = ",".join(uni[s].get(g, "") for g in genes)
        rows.append(r)
    os.makedirs(a.out, exist_ok=True)
    P.write_tsv(os.path.join(a.out, "orthogroups.tsv"), rows, fields)

    one = [g for g in groups if g["status"] == "1to1"]
    lines = [
        f"RBH edges: {edges:,}",
        f"components: {len(groups):,}  (1:1 {len(one):,}; conflict {len(groups) - len(one):,})",
        f"1:1 with all Plasmodium species: {sum(g['plasmodium_complete'] for g in one):,}",
        f"1:1 with all species: {sum(g['n_species'] == len(species) for g in one):,}",
    ]
    for s in species:
        lines.append(f"  {s}: in a 1:1 group {sum(s in g['members'] for g in one):,}")
    with open(os.path.join(a.out, "orthogroup_summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
