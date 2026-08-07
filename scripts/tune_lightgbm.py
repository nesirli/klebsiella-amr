#!/usr/bin/env python3
"""Optuna hyperparameter search for LightGBM resistance classifiers.

Searches on the training split using stratified k-fold CV and writes the best
hyperparameters (excluding data-dependent is_unbalance) as JSON. The training
scripts can then load this JSON via --params-input.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import optuna
import pandas as pd
from lightgbm import LGBMClassifier
from sklearn.model_selection import StratifiedKFold

optuna.logging.set_verbosity(optuna.logging.WARNING)


def gene_columns(df, all_antibiotics):
    non_gene = {"run", "collection_date", "year", "location", *all_antibiotics}
    return [c for c in df.columns if c not in non_gene]


def load_xy(csv_path, antibiotic, genes):
    df = pd.read_csv(csv_path).dropna(subset=[antibiotic])
    y = (df[antibiotic] == "R").astype(int).values
    x = df[genes].fillna(0).values.astype("float32")
    return x, y


def cv_score(x, y, hyperparams, n_splits=3, seed=42):
    if len(np.unique(y)) < 2:
        return 0.0
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    aucs = []
    for tr_idx, val_idx in skf.split(x, y):
        x_tr, x_val = x[tr_idx], x[val_idx]
        y_tr, y_val = y[tr_idx], y[val_idx]
        clf = LGBMClassifier(**{**hyperparams, "is_unbalance": True, "verbose": -1})
        clf.fit(x_tr, y_tr)
        proba = clf.predict_proba(x_val)[:, 1]
        if len(np.unique(y_val)) > 1:
            from sklearn.metrics import roc_auc_score
            aucs.append(roc_auc_score(y_val, proba))
    return float(np.mean(aucs)) if aucs else 0.0


def objective(trial, x, y):
    hyperparams = {
        "objective": "binary",
        "metric": "binary_logloss",
        "n_estimators": trial.suggest_int("n_estimators", 100, 500, log=True),
        "num_leaves": trial.suggest_int("num_leaves", 7, 63),
        "max_depth": trial.suggest_int("max_depth", 2, 8),
        "learning_rate": trial.suggest_float("learning_rate", 1e-3, 0.3, log=True),
        "subsample": trial.suggest_float("subsample", 0.5, 1.0),
        "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, 1.0),
        "min_child_samples": trial.suggest_int("min_child_samples", 2, 20),
        "reg_alpha": trial.suggest_float("reg_alpha", 1e-8, 10.0, log=True),
        "reg_lambda": trial.suggest_float("reg_lambda", 1e-8, 10.0, log=True),
        "random_state": 42,
    }
    return cv_score(x, y, hyperparams)


def parse_args():
    p = argparse.ArgumentParser(description="Tune LightGBM hyperparameters with Optuna")
    p.add_argument("--train-features", required=True)
    p.add_argument("--antibiotic", required=True)
    p.add_argument("--all-antibiotics", nargs="+", required=True)
    p.add_argument("--output", required=True, help="JSON file for best hyperparameters")
    p.add_argument("--n-trials", type=int, default=30)
    p.add_argument("--n-splits", type=int, default=3)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def main():
    args = parse_args()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    genes = gene_columns(pd.read_csv(args.train_features), args.all_antibiotics)
    x, y = load_xy(args.train_features, args.antibiotic, genes)

    if len(np.unique(y)) < 2 or len(y) < args.n_splits * 2:
        best = {"skipped": "not enough labeled data or class variation for tuning"}
        json.dump(best, open(args.output, "w"), indent=2)
        print(json.dumps(best, indent=2))
        return

    study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=args.seed))
    study.optimize(lambda trial: objective(trial, x, y), n_trials=args.n_trials, show_progress_bar=True)

    best = dict(study.best_params)
    best["n_trials"] = args.n_trials
    best["n_splits"] = args.n_splits
    best["best_cv_roc_auc"] = study.best_value
    json.dump(best, open(args.output, "w"), indent=2)
    print(json.dumps(best, indent=2))


if __name__ == "__main__":
    main()
