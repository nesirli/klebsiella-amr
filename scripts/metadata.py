#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import pandas as pd
import yaml


def parse_phenotype(ast_string, drug):
    if pd.isna(ast_string):
        return None
    match = re.search(rf"{re.escape(drug)}=([RS])", ast_string, re.IGNORECASE)
    return match.group(1).upper() if match else None


def load_config(path):
    c = yaml.safe_load(open(path))
    return {
        "input": c["metadata"],
        "delimiter": c.get("delimiter", ";"),
        "train_cutoff": c["splits"]["train_cutoff"],
        "test_years": c["splits"]["test_years"],
        "antibiotics": c["antibiotics"],
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse NCBI Pathogen Detection metadata into train/test splits"
    )
    parser.add_argument("--config", help="Path to config.yaml (alternative to explicit args)")
    parser.add_argument("--input")
    parser.add_argument("--delimiter", default=";")
    parser.add_argument("--train-output")
    parser.add_argument("--test-output")
    parser.add_argument("--samples-output")
    parser.add_argument("--train-cutoff", type=int)
    parser.add_argument("--test-years", type=int, nargs="+")
    parser.add_argument("--antibiotics", nargs="+")
    parser.add_argument("--max-samples", type=int, default=-1,
                        help="Cap each split at N samples; -1 means no cap")
    args = parser.parse_args()

    if args.config:
        cfg = load_config(args.config)
        for key, value in cfg.items():
            if getattr(args, key) is None:
                setattr(args, key, value)

    required = ["input", "train_output", "test_output", "samples_output",
                "train_cutoff", "test_years", "antibiotics"]
    for key in required:
        if getattr(args, key) is None:
            parser.error(f"--{key.replace('_', '-')} is required when --config is not provided")

    return args


def main():
    args = parse_args()

    Path(args.train_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.test_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.samples_output).parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep=args.delimiter, encoding="utf-8-sig")
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

    if args.max_samples > 0:
        train_df = train_df.head(args.max_samples)
        test_df = test_df.head(args.max_samples)

    train_df.to_csv(args.train_output, index=False)
    test_df.to_csv(args.test_output, index=False)

    samples = sorted(set(train_df["run"]).union(set(test_df["run"])))
    Path(args.samples_output).write_text("\n".join(samples) + "\n")

    mk_path = Path(args.samples_output).with_suffix(".mk")
    mk_path.write_text(f"SAMPLES := {' '.join(samples)}\n")

    print(f"train: {len(train_df)} samples | test: {len(test_df)} samples")


if __name__ == "__main__":
    main()
