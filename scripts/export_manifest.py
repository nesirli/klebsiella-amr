#!/usr/bin/env python3
"""Freeze everything the serving app needs to reproduce training conditions.

A model is only meaningful against the exact feature layout it was fitted on:
the same gene columns, in the same order, with the same decision threshold. All
three of those live in files `make report` deletes (the feature matrices) or
scatters (one metrics JSON per model per drug). This writes them to one file
under results/, which is kept.

Without it, a serving app has to guess the column order, and a silently
misaligned vector produces confident nonsense rather than an error.
"""
import argparse
import json
from datetime import date
from pathlib import Path

import pandas as pd

import common

SERVED_MODELS = ("xgboost", "lightgbm")
MODEL_FILE = {"xgboost": "{abx}_model.json", "lightgbm": "{abx}_model.txt"}


def parse_args():
    p = argparse.ArgumentParser(description="Write the serving manifest")
    p.add_argument("--train-features", required=True)
    p.add_argument("--models-dir", required=True)
    p.add_argument("--antibiotics", nargs="+", required=True)
    p.add_argument("--amrfinder-db", default="")
    p.add_argument("--output", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    models_dir = Path(args.models_dir)
    out_dir = Path(args.output).parent

    genes = common.gene_columns(pd.read_csv(args.train_features, nrows=0),
                               args.antibiotics)

    models = {}
    for name in SERVED_MODELS:
        per_abx = {}
        for abx in args.antibiotics:
            metrics_path = models_dir / name / f"{abx}_metrics.json"
            model_path = models_dir / name / MODEL_FILE[name].format(abx=abx)
            if not metrics_path.exists() or not model_path.exists():
                continue
            metrics = json.loads(metrics_path.read_text())
            if metrics.get("skipped"):
                continue
            # Relative to the manifest's own directory, so the same file works
            # from the repo (results/models/) and from the image (/app/artifacts/)
            # without rewriting. predict.py resolves against the manifest path.
            try:
                rel = model_path.resolve().relative_to(out_dir.resolve())
            except ValueError:
                rel = model_path
            per_abx[abx] = {
                # Fitted on the training split; serving with 0.5 instead would
                # discard it. Ranges from 0.31 to 0.73 across these models.
                "threshold": float(metrics.get("threshold", 0.5)),
                "model": str(rel),
                "test_roc_auc": metrics.get("roc_auc"),
                "test_balanced_accuracy": metrics.get("balanced_accuracy"),
            }
        if per_abx:
            models[name] = per_abx

    manifest = {
        "created": date.today().isoformat(),
        "genes": genes,
        "n_genes": len(genes),
        "antibiotics": list(args.antibiotics),
        "models": models,
        "amrfinder_db": args.amrfinder_db,
    }
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"manifest: {len(genes)} genes | "
          f"{sum(len(v) for v in models.values())} model/antibiotic pairs "
          f"-> {args.output}")


if __name__ == "__main__":
    main()
