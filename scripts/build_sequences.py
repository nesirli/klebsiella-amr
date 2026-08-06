#!/usr/bin/env python3
# /// script
# dependencies = ["pandas"]
# ///
"""Build a per-sample AMR-gene *nucleotide sequence* matrix for DNABERT-2.

The sequence-level analogue of build_features.py. Instead of a 0/1 presence
matrix, each cell holds the nucleotide sequence of that gene in that sample
(empty if absent). Genes are keyed by AMRFinder's Element symbol -- the SAME
key build_features.py uses -- so the DNABERT model sees exactly the genes the
tabular models do, but at sequence resolution (allelic variants/SNPs a 0/1
flag would collapse). Joined to train/test metadata like build_features.py.

AMRFinder --nucleotide_output FASTA header format (verified):
    >{contig}:{start}-{stop} strand:{+/-} {ELEMENT_SYMBOL} {description...}
so the 3rd whitespace field is the gene symbol.
"""
import argparse
from pathlib import Path

import pandas as pd


def genes_for_sample(fna_path):
    """Map gene symbol -> nucleotide sequence for one sample's FASTA.

    If a symbol occurs more than once (paralogs/duplicate hits), keep the
    longest sequence -- the most complete copy of the determinant.
    """
    seqs = {}
    symbol, chunks = None, []

    def flush():
        if symbol is None:
            return
        seq = "".join(chunks).upper()
        if len(seq) > len(seqs.get(symbol, "")):
            seqs[symbol] = seq

    with open(fna_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                parts = line[1:].split()
                # >contig:coords strand:+ SYMBOL description...
                symbol = parts[2] if len(parts) >= 3 else parts[-1]
                chunks = []
            else:
                chunks.append(line)
        flush()
    return seqs


def build_sequence_matrix(fna_files):
    accession_seqs = {}
    for f in fna_files:
        accession = Path(f).name.removesuffix("_amr_genes.fna")
        accession_seqs[accession] = genes_for_sample(f)

    all_genes = sorted({g for seqs in accession_seqs.values() for g in seqs})

    rows = []
    for accession, seqs in accession_seqs.items():
        row = {"run": accession}
        for gene in all_genes:
            row[gene] = seqs.get(gene, "")
        rows.append(row)

    return pd.DataFrame(rows), all_genes


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a per-sample AMR-gene nucleotide sequence matrix"
    )
    parser.add_argument("--seq-files", nargs="+", required=True)
    parser.add_argument("--train-metadata", required=True)
    parser.add_argument("--test-metadata", required=True)
    parser.add_argument("--train-output", required=True)
    parser.add_argument("--test-output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    Path(args.train_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.test_output).parent.mkdir(parents=True, exist_ok=True)

    seq_matrix, all_genes = build_sequence_matrix(args.seq_files)

    train_meta = pd.read_csv(args.train_metadata)
    test_meta = pd.read_csv(args.test_metadata)

    train_df = train_meta.merge(seq_matrix, on="run", how="inner")
    test_df = test_meta.merge(seq_matrix, on="run", how="inner")

    train_df.to_csv(args.train_output, index=False)
    test_df.to_csv(args.test_output, index=False)

    print(
        f"genes: {len(all_genes)} | "
        f"train: {len(train_df)} samples | test: {len(test_df)} samples"
    )


if __name__ == "__main__":
    main()
