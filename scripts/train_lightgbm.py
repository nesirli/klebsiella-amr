#!/usr/bin/env python3
"""Train + test + interpret one LightGBM resistance classifier per antibiotic.

Self-contained per-model script. Produces: saved model, hyperparameters,
test predictions, metrics, and SHAP feature importance (CSV + plot PNG).
"""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from lightgbm import LGBMClassifier
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
)

MODEL_NAME = "lightgbm"


def gene_columns(df, all_antibiotics):
    non_gene = {"run", "collection_date", "year", "location", *all_antibiotics}
    return [c for c in df.columns if c not in non_gene]


def load_xy(csv_path, antibiotic, genes):
    df = pd.read_csv(csv_path).dropna(subset=[antibiotic])
    y = (df[antibiotic] == "R").astype(int)
    x = df[genes].fillna(0)
    return x, y, df["run"]


def placeholder_png(path, text):
    fig, ax = plt.subplots(figsize=(6, 2))
    ax.text(0.5, 0.5, text, ha="center", va="center", wrap=True)
    ax.axis("off")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    p = argparse.ArgumentParser(description="Train a LightGBM resistance model")
    p.add_argument("--train-features", required=True)
    p.add_argument("--test-features", required=True)
    p.add_argument("--antibiotic", required=True)
    p.add_argument("--all-antibiotics", nargs="+", required=True)
    p.add_argument("--model-output", required=True)
    p.add_argument("--params-output", required=True)
    p.add_argument("--metrics-output", required=True)
    p.add_argument("--predictions-output", required=True)
    p.add_argument("--importance-output", required=True)
    p.add_argument("--shap-plot-output", required=True)
    p.add_argument("--params-input", help="Optional JSON with hyperparameters to override defaults")
    return p.parse_args()


def main():
    args = parse_args()
    for out in (args.model_output, args.params_output, args.metrics_output,
                args.predictions_output, args.importance_output, args.shap_plot_output):
        Path(out).parent.mkdir(parents=True, exist_ok=True)

    genes = gene_columns(pd.read_csv(args.train_features), args.all_antibiotics)
    x_train, y_train, _ = load_xy(args.train_features, args.antibiotic, genes)
    x_test, y_test, test_runs = load_xy(args.test_features, args.antibiotic, genes)

    base = {"model": MODEL_NAME, "antibiotic": args.antibiotic,
            "n_train": int(len(x_train)), "n_test": int(len(x_test)),
            "n_features": len(genes)}

    if len(x_train) == 0 or len(x_test) == 0 or y_train.nunique() < 2:
        base["skipped"] = ("need both R and S in train and >=1 test sample "
                           "(not enough labeled data for this antibiotic)")
        Path(args.model_output).write_bytes(b"")
        json.dump(base, open(args.params_output, "w"), indent=2)
        json.dump(base, open(args.metrics_output, "w"), indent=2)
        pd.DataFrame(columns=["run", "actual", "predicted", "probability_resistant"]
                     ).to_csv(args.predictions_output, index=False)
        pd.DataFrame(columns=["gene", "mean_abs_shap"]).to_csv(args.importance_output, index=False)
        placeholder_png(args.shap_plot_output, f"{args.antibiotic} ({MODEL_NAME}): skipped\n{base['skipped']}")
        print(json.dumps(base, indent=2))
        return

    n_pos = int(y_train.sum())
    n_neg = int(len(y_train) - n_pos)

    hyperparams = {
        "objective": "binary",
        "metric": "binary_logloss",
        "n_estimators": 200,
        "num_leaves": 15,
        "max_depth": 4,
        "learning_rate": 0.1,
        "subsample": 0.8,
        "colsample_bytree": 0.8,
        "min_child_samples": 2,
        "reg_alpha": 0.1,
        "reg_lambda": 0.1,
        # LightGBM's native imbalance handling: re-weights gradients so the
        # minority class contributes equally to training. More stable than
        # manual scale_pos_weight on tiny, skewed AMR splits.
        "is_unbalance": True,
        "random_state": 42,
        "verbose": -1,
    }
    if args.params_input:
        with open(args.params_input) as f:
            tuned = json.load(f)
        # Keep native imbalance handling unless the tuning run changed it.
        tuned.pop("is_unbalance", None)
        hyperparams.update(tuned)
    model = LGBMClassifier(**hyperparams)
    model.fit(x_train, y_train)
    model.booster_.save_model(args.model_output)

    params = {**base, "hyperparameters": hyperparams,
              "train_class_balance": {"R": n_pos, "S": n_neg}}
    json.dump(params, open(args.params_output, "w"), indent=2)

    y_pred = model.predict(x_test)
    y_proba = model.predict_proba(x_test)[:, 1]
    metrics = {**base,
               "accuracy": accuracy_score(y_test, y_pred),
               "balanced_accuracy": balanced_accuracy_score(y_test, y_pred),
               "f1": f1_score(y_test, y_pred, zero_division=0),
               "roc_auc": roc_auc_score(y_test, y_proba) if y_test.nunique() > 1 else None}
    json.dump(metrics, open(args.metrics_output, "w"), indent=2)
    pd.DataFrame({"run": test_runs, "actual": y_test.values,
                  "predicted": y_pred, "probability_resistant": y_proba}
                 ).to_csv(args.predictions_output, index=False)

    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(x_train)
    # LightGBM binary TreeExplainer may return a list [class0, class1]; take positive.
    if isinstance(shap_values, list):
        shap_values = shap_values[1]
    shap_values = np.asarray(shap_values)
    importance = (pd.DataFrame({"gene": genes,
                                "mean_abs_shap": np.abs(shap_values).mean(axis=0)})
                  .sort_values("mean_abs_shap", ascending=False).reset_index(drop=True))
    importance.to_csv(args.importance_output, index=False)

    shap.summary_plot(shap_values, x_train, show=False, max_display=20)
    plt.title(f"{args.antibiotic} ({MODEL_NAME}) SHAP")
    plt.savefig(args.shap_plot_output, dpi=150, bbox_inches="tight")
    plt.close()

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
