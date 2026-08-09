#!/usr/bin/env python3
"""Download paired-end reads from ENA.

ENA's filereport can return multiple FASTQ URLs per run (e.g. multiple lanes),
semicolon-separated. This script collects all R1 (_1.fastq) and R2 (_2.fastq)
URLs, downloads them, and concatenates them into the requested output files.
"""
import argparse
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


MAX_RETRIES = 5
RETRY_DELAY = 5


def _basename(url):
    return Path(url).name


def classify_fastq_urls(urls):
    """Split URLs into read1 / read2 based on standard ENA naming.

    Returns (r1_urls, r2_urls). Raises if the pairing is ambiguous or
    single-end only.
    """
    r1 = [u for u in urls if "_1.fastq" in _basename(u)]
    r2 = [u for u in urls if "_2.fastq" in _basename(u)]

    if not r1 and not r2 and len(urls) == 2:
        # Fallback: no standard naming, assume order is R1/R2.
        return [urls[0]], [urls[1]]

    if not r1 or not r2:
        raise ValueError(
            f"could not pair FASTQ URLs: {len(r1)} R1 and {len(r2)} R2 from {urls}"
        )

    return r1, r2


def fetch_fastq_urls(run_accession, retries=MAX_RETRIES, retry_delay=RETRY_DELAY):
    api_url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={run_accession}&result=read_run&fields=fastq_ftp&format=tsv"
    )
    req = urllib.request.Request(
        api_url,
        headers={"User-Agent": "klebsiella-amr-pipeline/1.0"},
    )
    last_error = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req) as resp:
                lines = resp.read().decode().strip().splitlines()
            if len(lines) < 2:
                raise ValueError("empty filereport")

            urls = []
            for line in lines[1:]:
                fields = line.split("\t")
                if len(fields) < 2 or not fields[1]:
                    continue
                urls.extend([u for u in fields[1].split(";") if u])

            if not urls:
                raise ValueError("no fastq_ftp URLs in filereport")

            r1, r2 = classify_fastq_urls(urls)
            return (
                [f"https://{u}" if not u.startswith("http") else u for u in r1],
                [f"https://{u}" if not u.startswith("http") else u for u in r2],
            )
        except Exception as e:
            last_error = e
            if attempt == retries:
                break
            print(
                f"retry {attempt}/{retries} for ENA filereport {run_accession}: {e}",
                file=sys.stderr,
            )
            time.sleep(retry_delay)
    raise last_error


def _remote_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req) as resp:
        return int(resp.headers.get("Content-Length", 0))


def download(url, dest, retries=MAX_RETRIES, retry_delay=RETRY_DELAY):
    """Resume-aware download of a single URL to dest."""
    attempt = 0
    while attempt < retries:
        attempt += 1
        try:
            existing = os.path.getsize(dest) if os.path.exists(dest) else 0
            remote = _remote_size(url)

            if existing >= remote > 0:
                print(f"{dest}: already complete ({existing} bytes)")
                return

            req = urllib.request.Request(url)
            if existing:
                req.add_header("Range", f"bytes={existing}-")

            with urllib.request.urlopen(req) as resp:
                mode = "ab" if existing and resp.status == 206 else "wb"
                if mode == "wb":
                    existing = 0
                with open(dest, mode) as f:
                    while chunk := resp.read(1024 * 1024):
                        f.write(chunk)

            actual = os.path.getsize(dest)
            if remote > 0 and actual != remote:
                raise IOError(f"incomplete download: got {actual} of {remote} bytes")
            return
        except urllib.error.HTTPError as e:
            if e.code == 416:
                # Local file is already complete or larger than remote; start over.
                print(f"{dest}: range not satisfiable, restarting download", file=sys.stderr)
                try:
                    os.remove(dest)
                except FileNotFoundError:
                    pass
                # Do not count this as a retry; the next loop iteration starts fresh.
                attempt -= 1
                time.sleep(retry_delay)
                continue
            if attempt == retries:
                raise
            print(f"retry {attempt}/{retries} for {url}: {e}", file=sys.stderr)
            time.sleep(retry_delay)
        except Exception as e:
            if attempt == retries:
                raise
            print(f"retry {attempt}/{retries} for {url}: {e}", file=sys.stderr)
            time.sleep(retry_delay)


def download_and_concat(urls, final_dest, retries=MAX_RETRIES, retry_delay=RETRY_DELAY):
    """Download one or more URLs and concatenate them into final_dest."""
    if len(urls) == 1:
        download(urls[0], final_dest, retries, retry_delay)
        return

    dest_dir = Path(final_dest).parent
    dest_dir.mkdir(parents=True, exist_ok=True)
    parts = []
    for i, url in enumerate(urls):
        part = dest_dir / f"{Path(final_dest).stem}.part{i}{Path(final_dest).suffix}"
        download(url, str(part), retries, retry_delay)
        parts.append(part)

    with open(final_dest, "wb") as out:
        for part in parts:
            with open(part, "rb") as f:
                while chunk := f.read(1024 * 1024):
                    out.write(chunk)

    for part in parts:
        part.unlink()


SKIP_MARKER_SUFFIX = ".skip"


def _try_fetch_fastq_urls(accession):
    try:
        return fetch_fastq_urls(accession)
    except urllib.error.HTTPError as e:
        if e.code in (403, 404):
            return None
        raise
    except Exception as e:
        msg = str(e).lower()
        if "could not pair fastq urls" in msg or "empty filereport" in msg:
            return None
        if "no fastq_ftp urls" in msg:
            return None
        raise


def main():
    parser = argparse.ArgumentParser(description="Download paired-end reads from ENA")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--out1", required=True)
    parser.add_argument("--out2", required=True)
    parser.add_argument(
        "--allow-skip",
        action="store_true",
        help="Create a .skip marker instead of failing for controlled-access or "
        "single-end samples",
    )
    args = parser.parse_args()

    result = _try_fetch_fastq_urls(args.accession)

    if result is None:
        if args.allow_skip:
            skip_marker = Path(args.out1).parent / f"{args.accession}{SKIP_MARKER_SUFFIX}"
            skip_marker.touch()
            print(
                f"SKIP: {args.accession} is not publicly accessible (controlled access, "
                f"missing data, or single-end only)",
                file=sys.stderr,
            )
            sys.exit(0)
        else:
            raise RuntimeError(
                f"{args.accession}: cannot fetch FASTQ URLs — sample may be "
                f"controlled-access or single-end. Use --allow-skip to skip gracefully."
            )

    r1_urls, r2_urls = result
    download_and_concat(r1_urls, args.out1)
    download_and_concat(r2_urls, args.out2)


if __name__ == "__main__":
    main()
