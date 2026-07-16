#!/usr/bin/env python3
import argparse
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
            urllib.request.urlretrieve(url, dest)
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
