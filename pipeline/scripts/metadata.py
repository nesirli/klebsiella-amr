import argparse
import re
from pathlib import Path

import pandas as pd


def parse_phenotype(ast_string, drug):
    if pd.isna(ast_string):
        return None
    match = re.search(rf"{re.escape(drug)}=([RS])", ast_string, re.IGNORECASE)
    return match.group(1).upper() if match else None


def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse NCBI Pathogen Detection metadata into train/test splits"
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--train-output", required=True)
    parser.add_argument("--test-output", required=True)
    parser.add_argument("--train-cutoff", type=int, required=True)
    parser.add_argument("--test-years", type=int, nargs="+", required=True)
    parser.add_argument("--antibiotics", nargs="+", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    Path(args.train_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.test_output).parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep=";", encoding="utf-8-sig")
    df.columns = (
        df.columns.str.strip()
        .str.replace("#", "", regex=False)
        .str.lower()
        .str.replace(" ", "_")
    )
    df["year"] = pd.to_numeric(df["collection_date"].astype(str).str[:4], errors="coerce")

    for drug in args.antibiotics:
        df[drug] = df["ast_phenotypes"].apply(lambda x, d=drug: parse_phenotype(x, d))

    df = df[["run", "collection_date", "year", "location", *args.antibiotics]]
    df = df.dropna(subset=["year"])

    train_df = df[df["year"] <= args.train_cutoff]
    test_df = df[df["year"].isin(args.test_years)]

    train_df.to_csv(args.train_output, index=False)
    test_df.to_csv(args.test_output, index=False)

    print(f"train: {len(train_df)} samples | test: {len(test_df)} samples")


if __name__ == "__main__":
    main()
