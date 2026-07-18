#!/usr/bin/env python3
# /// script
# dependencies = ["pandas"]
# ///
import argparse
from pathlib import Path

import pandas as pd


def genes_for_sample(amr_tsv):
    df = pd.read_csv(amr_tsv, sep="\t")
    return set(df["Element symbol"].dropna())


def build_gene_matrix(amr_files):
    accession_genes = {}
    for f in amr_files:
        accession = Path(f).name.removesuffix("_amr.tsv")
        accession_genes[accession] = genes_for_sample(f)

    all_genes = sorted(set.union(*accession_genes.values())) if accession_genes else []

    rows = []
    for accession, genes in accession_genes.items():
        row = {"run": accession}
        for gene in all_genes:
            row[gene] = int(gene in genes)
        rows.append(row)

    return pd.DataFrame(rows), all_genes


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a gene presence/absence feature matrix from amrfinder reports"
    )
    parser.add_argument("--amr-files", nargs="+", required=True)
    parser.add_argument("--train-metadata", required=True)
    parser.add_argument("--test-metadata", required=True)
    parser.add_argument("--train-output", required=True)
    parser.add_argument("--test-output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    Path(args.train_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.test_output).parent.mkdir(parents=True, exist_ok=True)

    gene_matrix, all_genes = build_gene_matrix(args.amr_files)

    train_meta = pd.read_csv(args.train_metadata)
    test_meta = pd.read_csv(args.test_metadata)

    train_df = train_meta.merge(gene_matrix, on="run", how="inner")
    test_df = test_meta.merge(gene_matrix, on="run", how="inner")

    train_df.to_csv(args.train_output, index=False)
    test_df.to_csv(args.test_output, index=False)

    print(
        f"genes: {len(all_genes)} | "
        f"train: {len(train_df)} samples | test: {len(test_df)} samples"
    )


if __name__ == "__main__":
    main()
