#!/usr/bin/env python3
"""Step 4 — align orthologs structurally and map every low-complexity
region onto its partners.

For each orthogroup with >= 2 models, all members are aligned to each other
with Foldseek (exhaustive within the group, 3Di+AA local alignment, so
inserted segments come out as gaps rather than being forced into register).
Then, for every fLPS region in every member A and every partner B:

  aligned_frac   fraction of the region's residues aligned to residues of B
  b_has_lcr      whether those B residues fall in one of B's own regions
  Placement uses only the CORE alignment: residue pairs with pLDDT >= 70 in
  both models. Disordered segments are placed relative to the ordered core
  around them, never by aligning disorder to disorder.
  location  core_shared         >= 50% of the region is core-aligned to B
            internal_insertion  region between core anchors; B's chain between
                                the anchors is short (<= 5 res or 25% of A's):
                                B lacks the segment
            internal_variable   between anchors, but B has its own segment there
            terminal_extension  outside the core, and B has no comparable tail
            terminal_variable   outside the core, B has a tail too
            poor_alignment      TM (shorter chain) < 0.3 or < 20 core pairs
  For internal regions, the residues flanking the region on each side that do
  align to B (anchors a1, a2 -> b1, b2) describe the insertion site in B:
  b_gap          residues of B between b1 and b2 (0 = B's chain is continuous
                 there: a clean insertion)
  b_anchor_dist  CA-CA distance b1-b2 in Angstrom
  b_site_sse     B's secondary structure at the site: 'loop' (both anchors in
                 coil), 'element' (both anchors inside the same helix/strand
                 run), 'boundary' (otherwise)
  Region descriptors from A's own model: mean pLDDT, SSE fractions.

Pairs whose overall structural alignment is poor (TM normalised by the
shorter chain < --min-tm) are written with location 'poor_alignment'.

Outputs: pairs.tsv, region_partner.tsv, regions.tsv (one row per region with
model-side descriptors), and _done/<og>.ok markers so the step resumes.
"""
import argparse
import concurrent.futures as cf
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import structlib as S  # noqa: E402

FMT = "query,target,qstart,qend,tstart,tend,qlen,tlen,alntmscore,qtmscore,ttmscore,lddt,prob,evalue,qaln,taln"
PAIR_F = ["og_id", "a_species", "a_gene", "b_species", "b_gene", "a_len", "b_len",
          "a_start", "a_end", "b_start", "b_end", "alntm", "qtm", "ttm", "tm_short",
          "lddt", "prob", "evalue", "aligned_res"]
RP_F = ["og_id", "a_species", "a_gene", "region_start", "region_end", "region_len",
        "lcr_class", "region_plddt", "b_species", "b_gene", "b_has_region_anywhere",
        "tm_short", "core_pairs", "aligned_frac", "b_has_lcr", "location", "b_gap",
        "b_gap_lcr_frac", "b_gap_plddt", "b_anchor_dist", "b_site_sse",
        "a_anchor_left", "a_anchor_right", "b_seg_start", "b_seg_end", "b_seg_plddt",
        "b_seg_class", "b_seg_class_overlap", "b_seg_n_frac", "b_seg_gs_frac", "b_seg_ked_frac"]
REG_F = ["og_id", "species", "gene", "uniprot", "region_start", "region_end", "region_len",
         "lcr_class", "signatures", "rel_position", "mean_plddt", "frac_plddt70",
         "frac_helix", "frac_strand", "protein_len", "protein_mean_plddt"]


def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def append_tsv(path, rows, fields, lock):
    with lock:
        new = not os.path.exists(path)
        with open(path, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                               extrasaction="ignore", lineterminator="\n")
            if new:
                w.writeheader()
            w.writerows(rows)


def run_foldseek(foldseek, files, tmp, threads):
    qdir = os.path.join(tmp, "q")
    os.makedirs(qdir)
    for f in files:
        os.symlink(os.path.abspath(f), os.path.join(qdir, os.path.basename(f)))
    out = os.path.join(tmp, "aln.m8")
    cmd = [foldseek, "easy-search", qdir, qdir, out, os.path.join(tmp, "t"),
           "--exhaustive-search", "1", "-e", "10", "--alignment-type", "2",
           "--format-output", FMT, "--threads", str(threads), "-v", "1"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(out):
        raise RuntimeError((r.stderr or r.stdout)[-400:])
    rows = []
    with open(out) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) == len(FMT.split(",")):
                rows.append(dict(zip(FMT.split(","), p)))
    return rows


def acc_of(name):
    # AF-<acc>-F1-model_v6(.pdb.gz) as Foldseek reports it
    base = os.path.basename(name)
    for suf in (".pdb.gz", ".pdb", ".gz"):
        if base.endswith(suf):
            base = base[: -len(suf)]
    return base.split("-")[1] if base.startswith("AF-") else base


def site_sse(sse, b1, b2):
    s1, s2 = sse[b1], sse[b2]
    if s1 == "C" and s2 == "C":
        return "loop"
    if s1 == s2 and s1 in "HE" and all(c == s1 for c in sse[b1:b2 + 1]):
        return "element"
    return "boundary"


def map_regions(regions, annot_seq, model_seq, min_identity, min_mapped=0.5):
    """Move fLPS region coordinates from the annotation sequence (VEuPathDB)
    onto the AlphaFold model sequence (UniProt). Returns (regions, identity,
    n_dropped). Regions whose residues mostly do not map are dropped; if the
    two sequences are too different, all regions are dropped."""
    if not regions:
        return [], (1.0 if annot_seq == model_seq else S.coord_map(annot_seq, model_seq)[1]), 0
    if not annot_seq:
        return [], 0.0, len(regions)
    cmap, ident = S.coord_map(annot_seq, model_seq)
    if ident < min_identity:
        return [], ident, len(regions)
    out, dropped = [], 0
    for r in regions:
        pos = [cmap[i] for i in range(r["start"] - 1, min(r["end"], len(cmap)))
               if cmap[i] is not None]
        if len(pos) < min_mapped * (r["end"] - r["start"] + 1):
            dropped += 1
            continue
        out.append({**r, "start": min(pos) + 1, "end": max(pos) + 1,
                    "annot_start": r["start"], "annot_end": r["end"]})
    return out, ident, dropped


def describe_b_segment(row, B, mb, j0, j1):
    """Fill the partner-segment columns for B residues j0..j1-1 (0-based)."""
    if j1 <= j0:
        return
    seq = mb["seq"][j0:j1]
    n = len(seq)
    row["b_seg_start"], row["b_seg_end"] = j0 + 1, j1
    row["b_seg_plddt"] = f"{sum(mb['plddt'][j0:j1]) / n:.1f}"
    row["b_seg_n_frac"] = f"{seq.count('N') / n:.3f}"
    row["b_seg_gs_frac"] = f"{(seq.count('G') + seq.count('S')) / n:.3f}"
    row["b_seg_ked_frac"] = f"{sum(seq.count(x) for x in 'KED') / n:.3f}"
    best, best_ov = "none", 0
    for r in B["regions"]:
        ov = min(j1, r["end"]) - max(j0, r["start"] - 1)
        if ov > best_ov:
            best, best_ov = r["lcr_class"], ov
    row["b_seg_class"] = best
    row["b_seg_class_overlap"] = best_ov


def process_og(og, members, a):
    """members: list of dicts with species, gene, uniprot, accession, path, regions."""
    models = {}
    mem_rows = []
    for m in members:
        mod = S.read_model(m["path"])
        mod["sse"] = S.ca_sse(mod["ca"])
        models[m["accession"]] = mod
        regs, ident, dropped = map_regions(m["regions"], a.annot.get((m["species"], m["gene"]), ""),
                                           mod["seq"], a.min_identity)
        m["regions"] = regs
        mem_rows.append({"og_id": og, "species": m["species"], "gene": m["gene"],
                         "accession": m["accession"], "model_len": len(mod["seq"]),
                         "annot_len": len(a.annot.get((m["species"], m["gene"]), "")),
                         "seq_identity": f"{ident:.3f}", "regions_mapped": len(regs),
                         "regions_dropped": dropped,
                         "mean_plddt": f"{sum(mod['plddt']) / max(len(mod['plddt']), 1):.1f}"})
    by_acc = {m["accession"]: m for m in members}
    tmp = tempfile.mkdtemp(prefix=f"{og}_", dir=a.tmp)
    try:
        hits = run_foldseek(a.foldseek, [m["path"] for m in members], tmp, a.threads_per_og)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # best alignment per ordered pair (query A, target B)
    best = {}
    for h in hits:
        qa, ta = acc_of(h["query"]), acc_of(h["target"])
        if qa == ta or qa not in by_acc or ta not in by_acc:
            continue
        key = (qa, ta)
        if key not in best or float(h["alntmscore"]) > float(best[key]["alntmscore"]):
            best[key] = h

    pair_rows, rp_rows, reg_rows = [], [], []
    for m in members:
        mod = models[m["accession"]]
        L = len(mod["seq"])
        for r in m["regions"]:
            s0, e0 = r["start"] - 1, r["end"]
            pl = mod["plddt"][s0:e0] or [0.0]
            ss = mod["sse"][s0:e0] or "C"
            reg_rows.append({
                "og_id": og, "species": m["species"], "gene": m["gene"],
                "uniprot": m["accession"], "region_start": r["start"],
                "region_end": r["end"], "region_len": e0 - s0,
                "lcr_class": r["lcr_class"], "signatures": r["signatures"],
                "rel_position": f"{((s0 + e0) / 2) / L:.3f}",
                "mean_plddt": f"{sum(pl) / len(pl):.1f}",
                "frac_plddt70": f"{sum(1 for x in pl if x >= 70) / len(pl):.3f}",
                "frac_helix": f"{ss.count('H') / len(ss):.3f}",
                "frac_strand": f"{ss.count('E') / len(ss):.3f}",
                "protein_len": L,
                "protein_mean_plddt": f"{sum(mod['plddt']) / L:.1f}",
            })

    for (qa, ta), h in best.items():
        A, B = by_acc[qa], by_acc[ta]
        ma, mb = models[qa], models[ta]
        pairs = list(S.aligned_pairs(h["qaln"], h["taln"], int(h["qstart"]), int(h["tstart"])))
        amap = {i: j for i, j in pairs}
        qtm, ttm = float(h["qtmscore"]), float(h["ttmscore"])
        tm_short = qtm if len(ma["seq"]) <= len(mb["seq"]) else ttm
        pair_rows.append({
            "og_id": og, "a_species": A["species"], "a_gene": A["gene"],
            "b_species": B["species"], "b_gene": B["gene"],
            "a_len": len(ma["seq"]), "b_len": len(mb["seq"]),
            "a_start": h["qstart"], "a_end": h["qend"], "b_start": h["tstart"],
            "b_end": h["tend"], "alntm": h["alntmscore"], "qtm": h["qtmscore"],
            "ttm": h["ttmscore"], "tm_short": f"{tm_short:.4f}", "lddt": h["lddt"],
            "prob": h["prob"], "evalue": h["evalue"], "aligned_res": len(pairs),
        })
        b_lcr = set()
        for r in B["regions"]:
            b_lcr.update(range(r["start"] - 1, r["end"]))
        # Core alignment: only residue pairs confidently modelled in BOTH
        # proteins. Disordered segments are never placed by aligning them to
        # each other (coil-to-coil alignment carries no homology signal);
        # they are placed relative to the ordered core that flanks them.
        pa, pb = ma["plddt"], mb["plddt"]
        core = {i: j for i, j in pairs if pa[i] >= a.core_plddt and pb[j] >= a.core_plddt}
        core_a = sorted(core)
        LA, LB = len(ma["seq"]), len(mb["seq"])
        for r in A["regions"]:
            s0, e0 = r["start"] - 1, r["end"]
            rl = e0 - s0
            mapped = [core[i] for i in range(s0, e0) if i in core]
            frac = len(mapped) / max(rl, 1)
            rpl = pa[s0:e0] or [0.0]
            row = {
                "og_id": og, "a_species": A["species"], "a_gene": A["gene"],
                "region_start": r["start"], "region_end": r["end"], "region_len": rl,
                "lcr_class": r["lcr_class"], "region_plddt": f"{sum(rpl) / len(rpl):.1f}",
                "b_species": B["species"], "b_gene": B["gene"],
                "b_has_region_anywhere": int(bool(B["regions"])),
                "tm_short": f"{tm_short:.4f}", "core_pairs": len(core),
                "aligned_frac": f"{frac:.3f}",
                "b_has_lcr": int(any(j in b_lcr for j in mapped)) if mapped else 0,
                "b_gap": "", "b_gap_lcr_frac": "", "b_gap_plddt": "",
                "b_anchor_dist": "", "b_site_sse": "",
                "a_anchor_left": "", "a_anchor_right": "",
                "b_seg_start": "", "b_seg_end": "", "b_seg_plddt": "", "b_seg_class": "",
                "b_seg_class_overlap": "", "b_seg_n_frac": "", "b_seg_gs_frac": "",
                "b_seg_ked_frac": "",
            }
            left = [i for i in core_a if i < s0]
            right = [i for i in core_a if i >= e0]
            if tm_short < a.min_tm or len(core) < a.min_core:
                row["location"] = "poor_alignment"
            elif frac >= 0.5:
                # region lies in the folded core and is aligned to B
                row["location"] = "core_shared"
            elif left and right:
                a1, a2 = left[-1], right[0]
                b1, b2 = core[a1], core[a2]
                row["a_anchor_left"], row["a_anchor_right"] = a1 + 1, a2 + 1
                if b2 <= b1:
                    row["location"] = "crossed"
                else:
                    gap = list(range(b1 + 1, b2))
                    g = len(gap)
                    row["b_gap"] = g
                    row["b_gap_lcr_frac"] = f"{sum(1 for j in gap if j in b_lcr) / g:.3f}" if g else "0"
                    row["b_gap_plddt"] = f"{sum(pb[j] for j in gap) / g:.1f}" if g else ""
                    row["b_anchor_dist"] = f"{S.dist(mb['ca'][b1], mb['ca'][b2]):.1f}"
                    row["b_site_sse"] = site_sse(mb["sse"], b1, b2)
                    describe_b_segment(row, B, mb, b1 + 1, b2)
                    # A's own span between anchors includes the region plus
                    # any linker; compare it with B's span
                    a_span = a2 - a1 - 1
                    if g <= max(a.ins_abs, a.ins_frac * a_span):
                        row["location"] = "internal_insertion"     # B lacks it
                    else:
                        row["location"] = "internal_variable"      # B has a segment there too
            else:
                # region outside the aligned core on one side: a terminal
                # extension. Compare A's tail with B's tail beyond the
                # nearest core anchor.
                if not core_a:
                    row["location"] = "poor_alignment"
                elif not left and not right:
                    # region spans the whole aligned core but is mostly not
                    # core-aligned itself: a long biased region containing
                    # the folded part
                    row["location"] = "core_overlap"
                else:
                    if not left:                                   # N-terminal side
                        anchor = right[0]
                        a_tail, b_tail = anchor, core[anchor]
                        describe_b_segment(row, B, mb, 0, core[anchor])
                    else:                                          # C-terminal side
                        anchor = left[-1]
                        a_tail, b_tail = LA - 1 - anchor, LB - 1 - core[anchor]
                        describe_b_segment(row, B, mb, core[anchor] + 1, LB)
                    row["b_gap"] = b_tail
                    row["a_anchor_left" if left else "a_anchor_right"] = anchor + 1
                    row["location"] = ("terminal_extension"
                                       if b_tail <= max(a.ins_abs, a.ins_frac * a_tail)
                                       else "terminal_variable")
            rp_rows.append(row)
    return pair_rows, rp_rows, reg_rows, mem_rows


MEM_F = ["og_id", "species", "gene", "accession", "model_len", "annot_len", "seq_identity",
         "regions_mapped", "regions_dropped", "mean_plddt"]


def load_annotation(seqs_dir):
    """(species, gene) -> primary-isoform annotation sequence."""
    from polylib import read_fasta, parse_veupath_header, primary_isoforms, is_pseudo
    out = {}
    for sp in ("pfalciparum", "pvivax", "pknowlesi", "pberghei", "tgondii"):
        p = os.path.join(seqs_dir, f"{sp}_proteins.fasta")
        if not os.path.exists(p):
            continue
        recs = []
        for h, seq in read_fasta(p):
            meta = parse_veupath_header(h)
            if not is_pseudo(meta):
                recs.append((meta, seq.upper().rstrip("*")))
        for g, (_, seq) in primary_isoforms(recs).items():
            out[(sp, g)] = seq
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", required=True, help="data/structures/_download_status.tsv")
    ap.add_argument("--structures", required=True, help="data/structures")
    ap.add_argument("--lcr-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--foldseek", default="foldseek")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--threads-per-og", type=int, default=1)
    ap.add_argument("--min-tm", type=float, default=0.3)
    ap.add_argument("--core-plddt", type=float, default=70.0,
                    help="both residues must reach this pLDDT to count as core-aligned")
    ap.add_argument("--min-core", type=int, default=20,
                    help="pairs with fewer core-aligned residues are poor alignments")
    ap.add_argument("--ins-abs", type=int, default=5,
                    help="B segment at most this long (or --ins-frac of A's) = insertion")
    ap.add_argument("--ins-frac", type=float, default=0.25)
    ap.add_argument("--only", default="", help="comma-separated og_ids (testing)")
    ap.add_argument("--seqs", required=True, help="data/seqs with <species>_proteins.fasta")
    ap.add_argument("--min-identity", type=float, default=0.9,
                    help="annotation vs model sequence identity needed to map regions")
    ap.add_argument("--tmp", default="")
    a = ap.parse_args()

    if not shutil.which(a.foldseek) and not os.path.exists(a.foldseek):
        sys.exit(f"ERROR: foldseek not found ({a.foldseek})")
    os.makedirs(os.path.join(a.out, "_done"), exist_ok=True)
    a.tmp = a.tmp or os.path.join(a.out, "_tmp")
    os.makedirs(a.tmp, exist_ok=True)

    a.annot = load_annotation(a.seqs)
    regions = defaultdict(list)
    for sp in ("pfalciparum", "pvivax", "pknowlesi", "pberghei", "tgondii"):
        p = os.path.join(a.lcr_dir, f"{sp}_lcr.tsv")
        if os.path.exists(p):
            for r in read_tsv(p):
                regions[(sp, r["gene"])].append({
                    "start": int(r["start"]), "end": int(r["end"]),
                    "lcr_class": r["lcr_class"], "signatures": r["signatures"]})

    ogs = defaultdict(list)
    for r in read_tsv(a.status):
        if not r.get("file"):
            continue
        ogs[r["og_id"]].append({
            "species": r["species"], "gene": r["gene"], "accession": r["accession"],
            "path": os.path.join(a.structures, r["file"]),
            "regions": regions.get((r["species"], r["gene"]), []),
        })
    todo = [og for og, m in sorted(ogs.items()) if len(m) >= 2]
    if a.only:
        keep = set(a.only.split(","))
        todo = [og for og in todo if og in keep]
    done = {f[:-3] for f in os.listdir(os.path.join(a.out, "_done"))}
    todo = [og for og in todo if og not in done]
    print(f"{len(todo):,} orthogroups to align ({len(done):,} already done)", flush=True)

    import threading
    lock = threading.Lock()
    paths = {k: os.path.join(a.out, f"{k}.tsv")
             for k in ("pairs", "region_partner", "regions", "members")}
    failed = 0

    def work(og):
        pr, rp, rg, mm = process_og(og, ogs[og], a)
        append_tsv(paths["pairs"], pr, PAIR_F, lock)
        append_tsv(paths["region_partner"], rp, RP_F, lock)
        append_tsv(paths["regions"], rg, REG_F, lock)
        append_tsv(paths["members"], mm, MEM_F, lock)
        open(os.path.join(a.out, "_done", og + ".ok"), "w").close()

    with cf.ThreadPoolExecutor(a.workers) as ex:
        futs = {ex.submit(work, og): og for og in todo}
        for i, f in enumerate(cf.as_completed(futs), 1):
            try:
                f.result()
            except Exception as exc:  # noqa: BLE001
                failed += 1
                print(f"  {futs[f]} failed: {exc}", flush=True)
            if i % 100 == 0 or i == len(todo):
                print(f"  {i:,}/{len(todo):,} done, {failed} failed", flush=True)
    shutil.rmtree(a.tmp, ignore_errors=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
