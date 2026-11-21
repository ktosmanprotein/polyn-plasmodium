#!/usr/bin/env python3
"""Step 1a — download annotated protein and CDS FASTA files from VEuPathDB.

One protein file and one CDS file per species, from the Current_Release
directory of PlasmoDB / ToxoDB. Resumable: files already on disk are kept.
The exact release file used is recorded in _sources.tsv so Methods can cite it.

If the automatic download fails (site layout changed, network blocked), the
script prints the page to download from by hand; drop the two files into the
sequence folder named <key>_proteins.fasta and <key>_cds.fasta and re-run.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polylib as P  # noqa: E402

UA = {"User-Agent": ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                     "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126 Safari/537.36 "
                     "polyn-pipeline/1.0"),
      "Accept": "*/*"}
LAST_ERRORS = []   # why probes failed, shown when nothing is found


def curl_ok():
    return shutil.which("curl") is not None


def curl_exists(url):
    r = subprocess.run(["curl", "-sS", "-L", "-o", "/dev/null", "-w", "%{http_code}",
                        "-r", "0-0", "-A", UA["User-Agent"], "--max-time", "60", url],
                       capture_output=True, text=True)
    code = r.stdout.strip()
    if code not in ("200", "206"):
        LAST_ERRORS.append(f"curl {code or 'error'} {r.stderr.strip()[:120]} <- {url}")
    return code in ("200", "206")


def http_get(url, timeout=120, retries=3):
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (attempt + 1))
    raise last


def find_file(listing_html, suffix):
    names = re.findall(r'href="([^"?/]+' + re.escape(suffix) + r')"', listing_html)
    names = sorted(set(names))
    return names[-1] if names else None


PROJECTS = {"plasmodb.org": "PlasmoDB", "toxodb.org": "ToxoDB",
            "cryptodb.org": "CryptoDB", "piroplasmadb.org": "PiroplasmaDB"}


def url_exists(url, timeout=30):
    """True if the server answers 200 to a HEAD (or a 1-byte GET) request."""
    for method, extra in (("HEAD", {}), ("GET", {"Range": "bytes=0-0"})):
        try:
            req = urllib.request.Request(url, headers={**UA, **extra}, method=method)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                if r.status in (200, 206):
                    return True
        except urllib.error.HTTPError as exc:
            LAST_ERRORS.append(f"python {method} HTTP {exc.code} <- {url}")
            if exc.code == 404:
                return False
        except Exception as exc:  # noqa: BLE001
            LAST_ERRORS.append(f"python {method} {type(exc).__name__}: {exc} <- {url}")
    # Python could not get an answer (TLS, proxy, bot filter): ask curl, which
    # uses the system's network stack and certificates.
    return curl_exists(url) if curl_ok() else False


def release_urls(host, org, suffix, project, releases):
    """Candidate file URLs in versioned release folders, newest first."""
    for n in releases:
        # Current layout (2026): files served by the raw-files service, as
        # linked from each site's Download data page. Older layout kept as a
        # fallback.
        for base in (f"https://{host}/a/service/raw-files/release-{n}/{org}/fasta/data/",
                     f"https://{host}/common/downloads/release-{n}/{org}/fasta/data/"):
            yield n, base, f"{project}-{n}_{org}{suffix}"


def locate(host, org, suffix, project, release_hint, exists=url_exists):
    """Find the newest release that has this file.

    Tries the release that already worked for this site first (all VEuPathDB
    sites share a release number), then counts down. Returns (url, filename,
    release) or None.
    """
    order = list(range(MAX_RELEASE, MIN_RELEASE - 1, -1))
    if release_hint in order:
        order.remove(release_hint)
        order.insert(0, release_hint)
    for n, base, fname in release_urls(host, org, suffix, project, order):
        if exists(base + fname):
            return base + fname, fname, n
    return None


MAX_RELEASE, MIN_RELEASE = 80, 60


def download(url, dest):
    tmp = dest + ".part"
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
            shutil.copyfileobj(r, fh, length=1 << 20)
    except Exception as exc:  # noqa: BLE001
        if not curl_ok():
            raise
        print(f"         python download failed ({exc}); retrying with curl")
        r = subprocess.run(["curl", "-sS", "-fL", "-A", UA["User-Agent"],
                            "--max-time", "1800", "-o", tmp, url],
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"curl failed: {r.stderr.strip()[:200]}")
    with open(tmp, "rb") as fh:
        head = fh.read(1)
    if head != b">":
        os.replace(tmp, dest + ".notfasta")
        raise RuntimeError(f"server returned something that is not FASTA; "
                           f"saved as {dest}.notfasta for inspection")
    os.replace(tmp, dest)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-free-mb", type=int, default=300)
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    free_mb = shutil.disk_usage(a.out).free // (1 << 20)
    if free_mb < a.min_free_mb:
        sys.exit(f"ERROR: only {free_mb} MB free at {a.out}; this step needs "
                 f"about 150 MB and refuses to run below {a.min_free_mb} MB.")
    print(f"Free disk at {a.out}: {free_mb:,} MB")

    sources_path = os.path.join(a.out, "_sources.tsv")
    sources = {r["key"] + r["kind"]: r for r in P.read_tsv(sources_path)} \
        if os.path.exists(sources_path) else {}

    failed = []
    release_hint = None
    for sp in P.read_tsv(a.species):
        key, host, org = sp["key"], sp["veupath_host"], sp["veupath_org"]
        project = sp.get("veupath_project") or PROJECTS.get(host, "")
        current = f"https://{host}/common/downloads/Current_Release/{org}/fasta/data/"
        listing = None
        for kind, suffix in (("proteins", "_AnnotatedProteins.fasta"),
                             ("cds", "_AnnotatedCDSs.fasta")):
            dest = os.path.join(a.out, f"{key}_{kind}.fasta")
            if os.path.exists(dest) and os.path.getsize(dest) > 1000:
                print(f"[cached] {key} {kind}")
                continue
            url = fname = None
            # 1. the Current_Release directory listing, where it still exists
            try:
                if listing is None:
                    listing = http_get(current, retries=1).decode("utf-8", "replace")
                fname = find_file(listing, suffix)
                if fname:
                    url = current + fname
            except Exception:  # noqa: BLE001
                listing = ""
            # 2. versioned release folders, newest first
            if not url:
                hit = locate(host, org, suffix, project, release_hint)
                if hit:
                    url, fname, release_hint = hit
            if not url:
                print(f"[FAIL]   {key} {kind}: not found under Current_Release "
                      f"or release-{MAX_RELEASE}..{MIN_RELEASE}")
                for e in [x for x in LAST_ERRORS if "release-71" in x][:4] or LAST_ERRORS[-4:]:
                    print(f"         {e}")
                LAST_ERRORS.clear()
                print(f"         download by hand from https://{host} (Download data)")
                print(f"         and save as {dest}")
                failed.append(f"{key}_{kind}")
                continue
            try:
                print(f"[get]    {key} {kind}: {fname}")
                download(url, dest)
                mb = os.path.getsize(dest) / (1 << 20)
                print(f"         {mb:.1f} MB")
                sources[key + kind] = {"key": key, "kind": kind, "file": fname,
                                       "url": url, "date": time.strftime("%Y-%m-%d")}
            except Exception as exc:  # noqa: BLE001
                print(f"[FAIL]   {key} {kind}: {exc} ({url})")
                failed.append(f"{key}_{kind}")

    P.write_tsv(sources_path, sorted(sources.values(), key=lambda r: (r["key"], r["kind"])),
                ["key", "kind", "file", "url", "date"])
    if failed:
        print(f"\n{len(failed)} file(s) missing: {', '.join(failed)}")
        return 1
    print("\nAll sequence files present.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
