#!/usr/bin/env python3
"""Aggregate per-model, per-antibiotic metrics into a single summary JSON."""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Summarize model results")
    parser.add_argument("--models-dir", required=True)
    parser.add_argument("--antibiotics", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    models = ["xgboost", "lightgbm", "nn"]
    rows = []
    for model in models:
        for abx in args.antibiotics:
            metrics_file = Path(args.models_dir) / model / f"{abx}_metrics.json"
            if not metrics_file.exists():
                continue
            data = json.load(open(metrics_file))
            rows.append({
                "model": model,
                "antibiotic": abx,
                "n_train": data.get("n_train"),
                "n_test": data.get("n_test"),
                "f1": data.get("f1"),
                "balanced_accuracy": data.get("balanced_accuracy"),
                "roc_auc": data.get("roc_auc"),
            })

    summary = {
        "models": models,
        "antibiotics": args.antibiotics,
        "metrics": rows,
    }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.output, "w"), indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
