#!/usr/bin/env python3
"""Download paired-end reads from ENA.

ENA's filereport can return multiple FASTQ URLs per run (e.g. multiple lanes),
semicolon-separated. This script collects all R1 (_1.fastq) and R2 (_2.fastq)
URLs, downloads them, and concatenates them into the requested output files.
"""
import argparse
import gzip
import os
import socket
import sys
import time
import urllib.error
import urllib.request
import zlib
from pathlib import Path


MAX_RETRIES = 5
RETRY_DELAY = 5

# Being unable to reach ENA says nothing about the accession, so a network
# outage must not spend the retry budget: the caller retires a sample whose
# download fails, and a brief outage once retired 57 good samples in a row.
# Network failures instead back off and wait, up to this many seconds.
NETWORK_OUTAGE_BUDGET = 1800
MAX_BACKOFF = 60

TRANSPORT_ERRORS = (socket.gaierror, ConnectionError, TimeoutError)
# ENA throttles with 403 as well as 429, and its gateways return 5xx under
# load; a 404 is the accession's own problem.
TRANSIENT_HTTP_CODES = (403, 429, 500, 502, 503, 504)


class PartialTransfer(IOError):
    """A truncated download that still advanced since the previous attempt.

    ENA closes long transfers early often enough that one large FASTQ can need
    several resumed passes. Each pass is progress, not failure, so it must not
    spend the retry budget that decides whether to retire the sample: a 186 MB
    file once climbed 77.7 -> 77.9 -> 78.0 MB across three passes and was
    retired mid-download for it.
    """


def is_transient(exc):
    """Whether retrying could plausibly succeed, vs. the accession being bad."""
    if isinstance(exc, PartialTransfer):
        return True
    if isinstance(exc, urllib.error.HTTPError):
        return exc.code in TRANSIENT_HTTP_CODES
    return (isinstance(exc, TRANSPORT_ERRORS)
            or isinstance(getattr(exc, "reason", None), TRANSPORT_ERRORS))


def with_retries(operation, describe, retries=MAX_RETRIES, retry_delay=RETRY_DELAY):
    """Call operation(), retrying on failure. -> whatever operation returns.

    Ordinary failures (a malformed filereport, a truncated body) get `retries`
    attempts. Network failures do not count against that budget at all: they
    back off exponentially and keep waiting until NETWORK_OUTAGE_BUDGET is
    spent, so a run that spans days survives the connection dropping.
    """
    attempt = 0
    waited = 0.0
    backoff = retry_delay
    while True:
        try:
            return operation()
        except PartialTransfer as e:
            # Forward progress. Reset the outage clock so a large file over a
            # flaky link keeps resuming for as long as it keeps advancing.
            print(f"{describe}: {e}; resuming", file=sys.stderr)
            waited = 0.0
            backoff = retry_delay
            time.sleep(retry_delay)
        except Exception as e:
            if is_transient(e):
                if waited >= NETWORK_OUTAGE_BUDGET:
                    raise
                pause = min(backoff, MAX_BACKOFF, NETWORK_OUTAGE_BUDGET - waited)
                print(f"cannot reach ENA ({e}); retrying {describe} in {pause:.0f}s "
                      f"[waited {waited:.0f}/{NETWORK_OUTAGE_BUDGET}s]", file=sys.stderr)
                time.sleep(pause)
                waited += pause
                backoff = min(backoff * 2, MAX_BACKOFF)
                continue
            attempt += 1
            if attempt >= retries:
                raise
            print(f"retry {attempt}/{retries} for {describe}: {e}", file=sys.stderr)
            time.sleep(retry_delay)


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
    def fetch():
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

    return with_retries(fetch, f"ENA filereport {run_accession}", retries, retry_delay)


def _gzip_ok(path):
    """Whether the file is a complete, readable gzip stream.

    Matching Content-Length is not proof the bytes are intact: one FASTQ
    arrived at exactly its advertised 98 MB and still ended mid-stream. Only
    fastp noticed, three steps later, by which point the sample was retired.
    """
    try:
        decompressed = 0
        with gzip.open(path, "rb") as fh:
            while chunk := fh.read(4 * 1024 * 1024):
                decompressed += len(chunk)
        # An empty file reads as a valid zero-member stream rather than
        # raising, so emptiness has to be rejected explicitly.
        return decompressed > 0
    except (OSError, EOFError, zlib.error):
        return False


def _remote_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req) as resp:
        return int(resp.headers.get("Content-Length", 0))


def download(url, dest, retries=MAX_RETRIES, retry_delay=RETRY_DELAY):
    """Resume-aware download of a single URL to dest."""
    def attempt():
        # Two passes at most: a 416 means the local file is already complete or
        # longer than the remote, and the retry after deleting it sends no
        # Range header, so it cannot 416 again.
        for _ in range(2):
            started_at = os.path.getsize(dest) if os.path.exists(dest) else 0
            existing = started_at
            remote = _remote_size(url)

            if existing >= remote > 0:
                print(f"{dest}: already complete ({existing} bytes)")
                return

            req = urllib.request.Request(url)
            if existing:
                req.add_header("Range", f"bytes={existing}-")

            try:
                with urllib.request.urlopen(req) as resp:
                    mode = "ab" if existing and resp.status == 206 else "wb"
                    if mode == "wb":
                        existing = 0
                    with open(dest, mode) as f:
                        while chunk := resp.read(1024 * 1024):
                            f.write(chunk)
            except urllib.error.HTTPError as e:
                if e.code == 416 and existing:
                    print(f"{dest}: range not satisfiable, restarting download",
                          file=sys.stderr)
                    try:
                        os.remove(dest)
                    except FileNotFoundError:
                        pass
                    continue
                raise

            actual = os.path.getsize(dest)
            if remote > 0 and actual != remote:
                if actual > started_at:
                    raise PartialTransfer(
                        f"got {actual} of {remote} bytes")
                raise IOError(
                    f"download stalled at {actual} of {remote} bytes")
            if not _gzip_ok(dest):
                # Right length, wrong bytes. Resuming would append to damage,
                # so discard and refetch from scratch. This counts against the
                # ordinary budget, so a file that is corrupt at the source is
                # retired rather than refetched forever.
                os.remove(dest)
                raise IOError(f"corrupt gzip despite a complete {actual} bytes")
            return
        raise IOError(f"{dest}: could not restart download after a range error")

    with_retries(attempt, url, retries, retry_delay)


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
