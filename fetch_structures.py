#!/usr/bin/env python3
"""Step 3 — download AlphaFold models for the structure set, stored gzipped.

Reads results/orthology/structure_set.tsv (from step 2c) and fetches one
AlphaFold DB model per UniProt accession into
  <out>/<species>/AF-<acc>-F1-model_v<N>.pdb.gz
Foldseek and the later steps read .pdb.gz directly; gzip cuts disk ~3x.

Resumable: a file already on disk is skipped. Models that other projects
already downloaded (--reuse folders) are copied, not fetched again. Tries model versions newest
first (6, 5, 4). Records every outcome in <out>/_download_status.tsv, so
Methods can state exactly which models were used and which were missing.
Proteins without a UniProt accession (mostly P. vivax Sal-1) are reported as
'no_uniprot' and skipped; they have no AlphaFold model.
"""
import argparse
import concurrent.futures as cf
import gzip
import os
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polylib as P  # noqa: E402

UA = {"User-Agent": "polyn-pipeline/1.0 (academic research)"}
VERSIONS = (6, 5, 4)


def url_for(acc, v):
    return f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{v}.pdb"


class RateLimited(Exception):
    """Server asked us to slow down (HTTP 429/503)."""
    def __init__(self, wait=0):
        super().__init__(f"rate limited (retry after {wait or '?'} s)")
        self.wait = wait


def get(url, timeout=60):
    """Return bytes, None for a 404, or raise for anything else."""
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read()
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        if exc.code in (429, 503):
            try:
                wait = float(exc.headers.get("Retry-After") or 0)
            except ValueError:
                wait = 0
            raise RateLimited(wait)
        raise
    except Exception:  # noqa: BLE001  (TLS/certificate trouble: try curl)
        if not shutil.which("curl"):
            raise
        r = subprocess.run(["curl", "-sS", "-L", "--max-time", str(timeout),
                            "-w", "\n%{http_code}", url], capture_output=True)
        body, _, code = r.stdout.rpartition(b"\n")
        if code.strip() == b"404":
            return None
        if code.strip() != b"200":
            raise RuntimeError(f"curl HTTP {code.decode().strip() or 'error'}")
        return body


def existing(d, acc):
    for v in VERSIONS:
        p = os.path.join(d, f"AF-{acc}-F1-model_v{v}.pdb.gz")
        if os.path.exists(p) and os.path.getsize(p) > 200:
            return p, v
    return None, None


LOCAL = {}   # accession -> (version, path) of a model already on disk elsewhere


def index_local(roots):
    """Find AlphaFold models already downloaded by other projects."""
    import re
    pat = re.compile(r"^AF-([A-Z0-9]+)-F1-model_v(\d+)\.pdb(\.gz)?$")
    for root in roots:
        if not os.path.isdir(root):
            print(f"  (reuse root not found, skipped: {root})")
            continue
        n = 0
        for dirpath, _, names in os.walk(root, followlinks=True):
            for nm in names:
                m = pat.match(nm)
                if not m:
                    continue
                p = os.path.join(dirpath, nm)
                v = int(m.group(2))
                if os.path.getsize(p) > 1000 and v >= LOCAL.get(m.group(1), (0, ""))[0]:
                    LOCAL[m.group(1)] = (v, p)
                    n += 1
        print(f"  indexed {n:,} models under {root}")


def copy_local(acc, d):
    v, src = LOCAL[acc]
    out = os.path.join(d, f"AF-{acc}-F1-model_v{v}.pdb.gz")
    if src.endswith(".gz"):
        shutil.copyfile(src, out + ".part")
    else:
        with open(src, "rb") as fi, gzip.open(out + ".part", "wb", compresslevel=6) as fo:
            shutil.copyfileobj(fi, fo)
    os.replace(out + ".part", out)
    return "reused", v, out


def fetch(acc, d, delay):
    path, v = existing(d, acc)
    if path:
        return "cached", v, path
    if acc in LOCAL:
        return copy_local(acc, d)
    last = "missing"
    for v in VERSIONS:
        for attempt in range(8):
            try:
                data = get(url_for(acc, v))
                break
            except RateLimited as exc:
                # back off hard: 10, 20, 40 ... s, or what the server asks
                last = f"error {exc}"
                time.sleep(max(exc.wait, 10 * 2 ** min(attempt, 5)))
                data = False
            except Exception as exc:  # noqa: BLE001
                last = f"error {exc}"
                time.sleep(2 * (attempt + 1))
                data = False
                if attempt >= 2:
                    break
        if data is False:
            return last, "", ""
        if data is None:
            continue
        if not data.startswith((b"HEADER", b"ATOM", b"REMARK", b"MODEL", b"CRYST1")):
            return "not_pdb", v, ""
        out = os.path.join(d, f"AF-{acc}-F1-model_v{v}.pdb.gz")
        with gzip.open(out + ".part", "wb", compresslevel=6) as fh:
            fh.write(data)
        os.replace(out + ".part", out)
        time.sleep(delay)
        return "downloaded", v, out
    return "not_in_afdb", "", ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, help="structure_set.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--delay", type=float, default=0.4)
    ap.add_argument("--min-free-gb", type=float, default=5.0)
    ap.add_argument("--limit", type=int, default=0, help="first N accessions only (test run)")
    ap.add_argument("--reuse", action="append", default=[],
                    help="folder holding AlphaFold models from other projects; "
                         "matching models are copied instead of downloaded (repeatable)")
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    free = shutil.disk_usage(a.out).free / 2**30
    print(f"Free disk at {a.out}: {free:,.1f} GB")
    if free < a.min_free_gb:
        sys.exit(f"ERROR: needs about 1.5 GB; refusing to start below {a.min_free_gb} GB free.")

    if a.reuse:
        print("Looking for models already on disk...")
        index_local(a.reuse)
    rows = P.read_tsv(a.set)
    jobs, status = [], []
    seen = set()
    for r in rows:
        acc = (r.get("uniprot") or "").split(",")[0].strip()
        if not acc:
            status.append({**r, "accession": "", "status": "no_uniprot", "version": "", "file": ""})
            continue
        if acc in seen:
            continue
        seen.add(acc)
        jobs.append((r, acc, os.path.join(a.out, r["species"])))
    if a.limit:
        jobs = jobs[:a.limit]
    for _, _, d in jobs:
        os.makedirs(d, exist_ok=True)
    reusable = sum(1 for _, acc, _ in jobs if acc in LOCAL)
    print(f"{len(jobs):,} models needed: {reusable:,} already on disk elsewhere, "
          f"{len(jobs) - reusable:,} to download; {len(status):,} proteins have no UniProt id")

    counts = {}
    t0 = time.time()
    with cf.ThreadPoolExecutor(a.workers) as ex:
        futs = {ex.submit(fetch, acc, d, a.delay): (r, acc) for r, acc, d in jobs}
        for i, f in enumerate(cf.as_completed(futs), 1):
            r, acc = futs[f]
            s, v, path = f.result()
            key = s.split()[0]
            counts[key] = counts.get(key, 0) + 1
            status.append({**r, "accession": acc, "status": s, "version": v,
                           "file": os.path.relpath(path, a.out) if path else ""})
            if i % 250 == 0 or i == len(jobs):
                rate = i / max(time.time() - t0, 1)
                eta = (len(jobs) - i) / rate / 60 if rate else 0
                print(f"  {i:,}/{len(jobs):,}  {counts}  ~{eta:.0f} min left", flush=True)

    fields = ["og_id", "species", "gene", "uniprot", "has_lcr", "accession",
              "status", "version", "file"]
    status.sort(key=lambda x: (x["og_id"], x["species"]))
    P.write_tsv(os.path.join(a.out, "_download_status.tsv"), status, fields)
    errors = sum(v for k, v in counts.items() if k in ("error", "not_pdb"))
    print(f"\nDone: {counts}")
    if errors:
        print(f"{errors} failed with errors; re-run to retry them (finished files are kept).")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
