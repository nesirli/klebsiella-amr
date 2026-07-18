#!/usr/bin/env python3
import argparse
import os
import sys
import time
import urllib.request


def fetch_fastq_urls(run_accession):
    api_url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={run_accession}&result=read_run&fields=fastq_ftp&format=tsv"
    )
    with urllib.request.urlopen(api_url) as resp:
        lines = resp.read().decode().strip().splitlines()
    fastq_ftp = lines[1].split("\t")[1]
    url1, url2 = fastq_ftp.split(";")
    return f"https://{url1}", f"https://{url2}"


def download(url, dest, retries=5, retry_delay=5):
    for attempt in range(1, retries + 1):
        try:
            existing = os.path.getsize(dest) if os.path.exists(dest) else 0
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

            expected = existing + int(resp.headers.get("Content-Length", 0))
            actual = os.path.getsize(dest)
            if expected and actual != expected:
                raise IOError(f"incomplete download: got {actual} of {expected} bytes")
            return
        except Exception as e:
            if attempt == retries:
                raise
            print(f"retry {attempt}/{retries} for {url}: {e}", file=sys.stderr)
            time.sleep(retry_delay)


def main():
    parser = argparse.ArgumentParser(description="Download paired-end reads from ENA")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--out1", required=True)
    parser.add_argument("--out2", required=True)
    args = parser.parse_args()

    url1, url2 = fetch_fastq_urls(args.accession)
    download(url1, args.out1)
    download(url2, args.out2)


if __name__ == "__main__":
    main()
