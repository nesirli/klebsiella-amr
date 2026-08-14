#!/usr/bin/env python3
"""Shared plumbing for the per-model train and tune scripts.

Every model follows the same contract: read the gene matrix built by
build_features.py (or build_sequences.py), fit one binary R-vs-S classifier for
one antibiotic, and emit the same six artifacts (model, params, metrics,
predictions, importance CSV, importance plot). Only the estimator and its
importance method differ between models; everything around them lives here.

Tuner/trainer handoff
---------------------
Tuners write {"hyperparameters": {...}, "tuning": {...}} so that search
bookkeeping (n_trials, best CV score) can never be mistaken for a model
hyperparameter. Trainers read only the "hyperparameters" block. The full dict
is rebuilt by replaying the tuner's own suggest() through an Optuna FixedTrial,
so a structured hyperparameter (the MLP's `hidden` layer list, assembled from
two scalar suggestions) survives the round trip intact.
"""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)

# Columns in a features CSV that are metadata rather than gene features. The
# antibiotic label columns are excluded separately since they are configurable.
METADATA_COLUMNS = ("run", "collection_date", "year", "location")

PREDICTION_COLUMNS = ["run", "actual", "predicted", "probability_resistant"]
IMPORTANCE_COLUMNS = ["gene", "importance"]

SKIP_REASON = ("need both R and S in train and >=1 test sample "
               "(not enough labeled data for this antibiotic)")
TUNE_SKIP_REASON = "not enough labeled data or class variation for tuning"


# Data ------------------------------------------------------------------------

def gene_columns(df, all_antibiotics):
    """Feature column names: everything that is neither metadata nor a label."""
    non_gene = {*METADATA_COLUMNS, *all_antibiotics}
    return [c for c in df.columns if c not in non_gene]

def read_genes(features_csv, all_antibiotics):
    """Gene columns of a features CSV, read from the header alone."""
    return gene_columns(pd.read_csv(features_csv, nrows=0), all_antibiotics)


def load_labelled(features_csv, antibiotic, dtype=None):
    """Rows of a features CSV that carry an R/S call for this antibiotic."""
    df = pd.read_csv(features_csv, dtype=dtype).dropna(subset=[antibiotic])
    return df.reset_index(drop=True)


def load_xy(features_csv, antibiotic, genes):
    """-> (X frame, y 0/1 series, run series) for the labelled rows."""
    df = load_labelled(features_csv, antibiotic)
    y = (df[antibiotic] == "R").astype(int)
    return df[genes].fillna(0), y, df["run"]


def load_arrays(features_csv, antibiotic, genes):
    """load_xy as numpy, for the torch models and the tuners."""
    x, y, runs = load_xy(features_csv, antibiotic, genes)
    return x.to_numpy(dtype="float32"), y.to_numpy(), runs


def class_weights(y):
    """Balanced CrossEntropyLoss weights [w_susceptible, w_resistant]."""
    import torch

    n_pos = int(np.sum(y))
    n_neg = int(len(y) - n_pos)
    return torch.tensor([len(y) / (2 * max(n_neg, 1)),
                         len(y) / (2 * max(n_pos, 1))], dtype=torch.float32)


def enough_to_fit(y_train, n_test, min_train=1):
    """Whether there is both an R and an S to learn from and something to score."""
    return (len(y_train) >= min_train and n_test > 0
            and len(np.unique(np.asarray(y_train))) > 1)


# Artifacts -------------------------------------------------------------------

def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2) + "\n")


def print_json(obj):
    """Echo an artifact to stdout so it lands in the make log."""
    print(json.dumps(obj, indent=2))


def write_predictions(path, runs, y_true, y_pred, y_proba):
    pd.DataFrame({"run": np.asarray(runs), "actual": np.asarray(y_true),
                  "predicted": np.asarray(y_pred),
                  "probability_resistant": np.asarray(y_proba)}
                 ).to_csv(path, index=False)


def write_importance(path, genes, values):
    """Gene importance CSV, most important first. -> the written frame."""
    df = (pd.DataFrame({"gene": list(genes), "importance": list(values)})
          .sort_values("importance", ascending=False).reset_index(drop=True))
    df.to_csv(path, index=False)
    return df


def importance_plot(path, importance, title, xlabel, top_n=20):
    top = importance.head(top_n)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(top["gene"][::-1], top["importance"][::-1])
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def placeholder_png(path, text):
    fig, ax = plt.subplots(figsize=(6, 2))
    ax.text(0.5, 0.5, text, ha="center", va="center", wrap=True)
    ax.axis("off")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def compute_metrics(base, y_true, y_pred, y_proba):
    """The metric block every model reports, with `base` merged in.

    roc_auc and pr_auc are undefined for a single-class test split, and are
    reported as null rather than omitted so the schema stays stable.
    """
    y_true = np.asarray(y_true)
    both_classes = len(np.unique(y_true)) > 1
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return {
        **base,
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
        "roc_auc": float(roc_auc_score(y_true, y_proba)) if both_classes else None,
        "pr_auc": float(average_precision_score(y_true, y_proba)) if both_classes else None,
        "confusion_matrix": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)},
    }


def write_skipped(args, base, model_name, reason=SKIP_REASON):
    """Emit well-formed empty artifacts for an antibiotic we cannot model.

    The pipeline runs one target per (model, antibiotic) pair; a pair with too
    little labelled data still has to produce every output file or make would
    treat the whole run as failed.
    """
    base = {**base, "skipped": reason}
    Path(args.model_output).write_bytes(b"")
    write_json(args.params_output, base)
    write_json(args.metrics_output, base)
    pd.DataFrame(columns=PREDICTION_COLUMNS).to_csv(args.predictions_output, index=False)
    pd.DataFrame(columns=IMPORTANCE_COLUMNS).to_csv(args.importance_output, index=False)
    placeholder_png(args.importance_plot_output,
                    f"{args.antibiotic} ({model_name}): skipped\n{reason}")
    print_json(base)


# Command line ----------------------------------------------------------------

def train_parser(description):
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--train-features", required=True)
    p.add_argument("--test-features", required=True)
    p.add_argument("--antibiotic", required=True)
    p.add_argument("--all-antibiotics", nargs="+", required=True)
    p.add_argument("--model-output", required=True)
    p.add_argument("--params-output", required=True)
    p.add_argument("--metrics-output", required=True)
    p.add_argument("--predictions-output", required=True)
    p.add_argument("--importance-output", required=True)
    p.add_argument("--importance-plot-output", required=True)
    p.add_argument("--params-input",
                   help="Tuning JSON from the matching tune_*.py; its "
                        "'hyperparameters' block overrides the defaults")
    return p


def prepare_outputs(args):
    """Create the parent directory of every --*-output path."""
    for name, value in vars(args).items():
        if name.endswith("_output") and value:
            Path(value).parent.mkdir(parents=True, exist_ok=True)


# Tuning ----------------------------------------------------------------------

def tune_parser(description, default_trials):
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--train-features", required=True)
    p.add_argument("--antibiotic", required=True)
    p.add_argument("--all-antibiotics", nargs="+", required=True)
    p.add_argument("--output", required=True, help="JSON file for best hyperparameters")
    p.add_argument("--n-trials", type=int, default=default_trials)
    p.add_argument("--n-splits", type=int, default=3)
    p.add_argument("--seed", type=int, default=42)
    return p


def write_tuning_result(path, hyperparameters, tuning):
    """Write the tuner -> trainer handoff file.

    Hyperparameters and search bookkeeping are kept in separate blocks so that
    the trainer can splat the former into an estimator constructor without
    also passing it n_trials or the best CV score.
    """
    doc = {"hyperparameters": hyperparameters, "tuning": tuning}
    write_json(path, doc)
    print_json(doc)


def load_tuned(params_input):
    """Hyperparameters block of a tuning JSON; {} if absent or skipped."""
    if not params_input:
        return {}
    with open(params_input) as f:
        return json.load(f).get("hyperparameters") or {}


def cv_roc_auc(x, y, fit_predict, n_splits, seed=42):
    """Mean ROC AUC over stratified folds.

    `fit_predict(x_train, y_train, x_val)` returns P(resistant) for x_val.
    Folds whose validation split ends up single-class contribute no AUC.
    """
    from sklearn.model_selection import StratifiedKFold

    if len(np.unique(y)) < 2:
        return 0.0
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    aucs = []
    for train_idx, val_idx in skf.split(x, y):
        y_val = y[val_idx]
        if len(np.unique(y_val)) < 2:
            continue
        proba = fit_predict(x[train_idx], y[train_idx], x[val_idx])
        aucs.append(roc_auc_score(y_val, proba))
    return float(np.mean(aucs)) if aucs else 0.0


def run_tuning(args, suggest, fit_predict):
    """Shared Optuna driver.

    `suggest(trial)` builds a complete hyperparameter dict from a trial;
    `fit_predict(hyperparams, x_train, y_train, x_val)` returns validation
    probabilities. The winning dict is rebuilt by replaying `suggest` over the
    best trial's parameters, so whatever shape `suggest` produces -- including
    values assembled from several suggestions -- is what gets written out.
    """
    import optuna

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    genes = read_genes(args.train_features, args.all_antibiotics)
    x, y, _ = load_arrays(args.train_features, args.antibiotic, genes)

    if len(np.unique(y)) < 2 or len(y) < args.n_splits * 2:
        write_tuning_result(args.output, {}, {"skipped": TUNE_SKIP_REASON})
        return

    def objective(trial):
        hyperparams = suggest(trial)
        return cv_roc_auc(
            x, y,
            lambda x_tr, y_tr, x_val: fit_predict(hyperparams, x_tr, y_tr, x_val),
            args.n_splits, args.seed,
        )

    study = optuna.create_study(direction="maximize",
                                sampler=optuna.samplers.TPESampler(seed=args.seed))
    study.optimize(objective, n_trials=args.n_trials, show_progress_bar=True)

    write_tuning_result(
        args.output,
        suggest(optuna.trial.FixedTrial(study.best_params)),
        {"n_trials": args.n_trials, "n_splits": args.n_splits,
         "best_cv_roc_auc": float(study.best_value)},
    )
