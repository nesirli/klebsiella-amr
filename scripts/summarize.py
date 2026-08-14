#!/usr/bin/env python3
"""Aggregate per-model, per-antibiotic metrics and top features into one JSON."""
import argparse
import json
from pathlib import Path

import pandas as pd

# dnabert is optional (`make dnabert`); models with no metrics on disk are
# simply left out of the summary.
MODELS = ["xgboost", "lightgbm", "nn", "dnabert"]

METRIC_KEYS = ["n_train", "n_test", "accuracy", "balanced_accuracy", "f1",
               "precision", "recall", "roc_auc", "pr_auc", "confusion_matrix",
               "skipped"]


def top_features(importance_csv, top_n):
    """The highest-importance genes for one model/antibiotic, if recorded."""
    if not importance_csv.exists():
        return []
    df = pd.read_csv(importance_csv)
    if df.empty:
        return []
    return [{"gene": row.gene, "importance": float(row.importance)}
            for row in df.head(top_n).itertuples()]


def main():
    parser = argparse.ArgumentParser(description="Summarize model results")
    parser.add_argument("--models-dir", required=True)
    parser.add_argument("--antibiotics", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--top-n", type=int, default=10,
                        help="Genes to carry over per model/antibiotic")
    args = parser.parse_args()

    models_dir = Path(args.models_dir)
    rows, present = [], []
    for model in MODELS:
        metrics_files = {abx: models_dir / model / f"{abx}_metrics.json"
                         for abx in args.antibiotics}
        if not any(f.exists() for f in metrics_files.values()):
            continue
        present.append(model)
        for abx, metrics_file in metrics_files.items():
            if not metrics_file.exists():
                continue
            data = json.loads(metrics_file.read_text())
            rows.append({
                "model": model,
                "antibiotic": abx,
                **{key: data.get(key) for key in METRIC_KEYS},
                "top_features": top_features(
                    models_dir / model / f"{abx}_importance.csv", args.top_n),
            })

    summary = {"models": present, "antibiotics": args.antibiotics, "metrics": rows}

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
