#!/usr/bin/env python3
"""Compute minimal assembly statistics (N50, contigs, GC%, total length)."""
import argparse
from pathlib import Path


def parse_fasta(path):
    seqs = []
    current = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    seqs.append("".join(current))
                    current = []
            else:
                current.append(line.upper())
        if current:
            seqs.append("".join(current))
    return seqs


def n50(lengths):
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    half = total / 2
    cum = 0
    for length in lengths:
        cum += length
        if cum >= half:
            return length
    return 0


def main():
    parser = argparse.ArgumentParser(description="Compute assembly statistics")
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    seqs = parse_fasta(args.assembly)
    lengths = [len(s) for s in seqs]
    total = sum(lengths)
    gc = sum(s.count("G") + s.count("C") for s in seqs)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as fh:
        fh.write("sample\tcontigs\ttotal_length\tn50\tgc_percent\n")
        fh.write(f"{Path(args.assembly).stem}\t{len(seqs)}\t{total}\t{n50(lengths)}\t"
                 f"{100.0 * gc / total if total else 0.0:.2f}\n")


if __name__ == "__main__":
    main()
