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


def stratified_sample(df, n, antibiotics, seed=42):
    """Take n rows, preserving the joint R/S profile across all antibiotics.

    A plain .head(n) takes whatever rows the input happens to list first, and
    NCBI metadata arrives grouped by submitting study: capping at 20 once drew
    18 R and 1 S for ceftazidime out of a pool that is 51% R, which left the
    tuners nothing to learn from. Rows are grouped by their full resistance
    pattern ("R|S|NA|R"), each group contributes its proportional share
    (largest remainder for the leftovers), and the draw within a group is
    shuffled under a fixed seed so the choice is representative but repeatable.
    """
    if n <= 0 or n >= len(df):
        return df

    shuffled = df.sample(frac=1, random_state=seed)
    pattern = shuffled[antibiotics].fillna("NA").agg("|".join, axis=1)

    quota = pattern.value_counts(normalize=True) * n
    take = quota.astype(int)
    shortfall = n - take.sum()
    if shortfall:
        leftovers = (quota - take).sort_values(ascending=False)
        take[leftovers.index[:shortfall]] += 1

    picked = [shuffled[pattern == key].head(count)
              for key, count in take.items() if count]
    return pd.concat(picked).sort_index()


def random_split(df, antibiotics, test_fraction, seed=42):
    """Stratified random train/test split, as a control for the temporal one.

    The temporal split is confounded with which studies deposited isolates
    when: ceftazidime is 47% resistant before 2021 and 85% after, and every
    model scores near chance on it (test ROC AUC 0.53-0.66) despite reaching
    0.94 in cross-validation within the training years. Splitting at random,
    stratified on the joint resistance profile, keeps both sides on the same
    distribution and so measures how much signal the features carry, separately
    from how well it survives a shift in time.
    """
    pattern = df[antibiotics].fillna("NA").agg("|".join, axis=1)
    shuffled = df.sample(frac=1, random_state=seed)
    pattern = pattern.loc[shuffled.index]

    test_idx = []
    for key in pattern.unique():
        group = shuffled.index[pattern == key]
        n_test = int(round(len(group) * test_fraction))
        # Keep singleton profiles in train rather than spending them on test.
        test_idx.extend(group[:n_test])

    test_set = set(test_idx)
    is_test = shuffled.index.isin(test_set)
    return shuffled[~is_test].sort_index(), shuffled[is_test].sort_index()


def load_config(path):
    c = yaml.safe_load(open(path))
    return {
        "input": c["metadata"],
        "delimiter": c.get("delimiter", ";"),
        "train_cutoff": c["splits"]["train_cutoff"],
        "test_years": c["splits"]["test_years"],
        "split_mode": c["splits"].get("mode", "temporal"),
        "test_fraction": c["splits"].get("test_fraction", 0.2),
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
    parser.add_argument("--split-mode", choices=["temporal", "random"],
                        help="temporal: train on <=cutoff, test on test_years. "
                             "random: stratified random split (control for "
                             "distribution shift between the two eras)")
    parser.add_argument("--test-fraction", type=float,
                        help="Test share when --split-mode=random")
    args = parser.parse_args()

    if args.config:
        cfg = load_config(args.config)
        for key, value in cfg.items():
            if getattr(args, key) is None:
                setattr(args, key, value)

    args.split_mode = args.split_mode or "temporal"
    args.test_fraction = args.test_fraction if args.test_fraction is not None else 0.2

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

    # An isolate with no R/S call for any modelled antibiotic is dropped from
    # every model at fit time, so downloading and assembling it buys nothing.
    labelled = df[args.antibiotics].notna().any(axis=1)
    if not labelled.all():
        print(f"dropping {(~labelled).sum()} samples with no R/S label "
              f"for any of: {', '.join(args.antibiotics)}")
        df = df[labelled]

    if args.split_mode == "random":
        train_df, test_df = random_split(df, args.antibiotics, args.test_fraction)
    else:
        train_df = df[df["year"] <= args.train_cutoff]
        test_df = df[df["year"].isin(args.test_years)]

    if args.max_samples > 0:
        train_df = stratified_sample(train_df, args.max_samples, args.antibiotics)
        test_df = stratified_sample(test_df, args.max_samples, args.antibiotics)

    train_df.to_csv(args.train_output, index=False)
    test_df.to_csv(args.test_output, index=False)

    samples = sorted(set(train_df["run"]).union(set(test_df["run"])))
    Path(args.samples_output).write_text("\n".join(samples) + "\n")

    mk_path = Path(args.samples_output).with_suffix(".mk")
    mk_path.write_text(f"SAMPLES := {' '.join(samples)}\n")

    print(f"{args.split_mode} split | train: {len(train_df)} samples "
          f"| test: {len(test_df)} samples")


if __name__ == "__main__":
    main()
