#!/usr/bin/env python3
# /// script
# dependencies = ["pandas", "scikit-learn", "xgboost"]
# ///
import argparse
import json
from pathlib import Path

import pandas as pd
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from xgboost import XGBClassifier


def load_xy(csv_path, antibiotic, gene_columns):
    df = pd.read_csv(csv_path)
    df = df.dropna(subset=[antibiotic])
    y = (df[antibiotic] == "R").astype(int)
    x = df[gene_columns].fillna(0)
    return x, y, df["run"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train an XGBoost resistance classifier for one antibiotic"
    )
    parser.add_argument("--train-features", required=True)
    parser.add_argument("--test-features", required=True)
    parser.add_argument("--antibiotic", required=True)
    parser.add_argument("--all-antibiotics", nargs="+", required=True)
    parser.add_argument("--metrics-output", required=True)
    parser.add_argument("--predictions-output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    non_gene_columns = {
        "run", "collection_date", "year", "location", *args.all_antibiotics,
    }
    train_df = pd.read_csv(args.train_features)
    gene_columns = [c for c in train_df.columns if c not in non_gene_columns]

    x_train, y_train, _ = load_xy(args.train_features, args.antibiotic, gene_columns)
    x_test, y_test, test_runs = load_xy(args.test_features, args.antibiotic, gene_columns)

    Path(args.metrics_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.predictions_output).parent.mkdir(parents=True, exist_ok=True)

    if len(x_train) == 0 or len(x_test) == 0 or y_train.nunique() < 2:
        metrics = {
            "antibiotic": args.antibiotic,
            "n_train": len(x_train),
            "n_test": len(x_test),
            "skipped": "not enough labeled train/test data for this antibiotic "
                       "(need >=1 sample of each class in train, >=1 test sample)",
        }
        with open(args.metrics_output, "w") as f:
            json.dump(metrics, f, indent=2)
        pd.DataFrame(columns=["run", "actual", "predicted", "probability_resistant"]).to_csv(
            args.predictions_output, index=False
        )
        print(json.dumps(metrics, indent=2))
        return

    model = XGBClassifier(eval_metric="logloss")
    model.fit(x_train, y_train)

    y_pred = model.predict(x_test)
    y_proba = model.predict_proba(x_test)[:, 1]

    metrics = {
        "antibiotic": args.antibiotic,
        "n_train": len(x_train),
        "n_test": len(x_test),
        "accuracy": accuracy_score(y_test, y_pred),
        "f1": f1_score(y_test, y_pred),
        "roc_auc": roc_auc_score(y_test, y_proba) if y_test.nunique() > 1 else None,
    }

    with open(args.metrics_output, "w") as f:
        json.dump(metrics, f, indent=2)

    predictions = pd.DataFrame({
        "run": test_runs,
        "actual": y_test.values,
        "predicted": y_pred,
        "probability_resistant": y_proba,
    })
    predictions.to_csv(args.predictions_output, index=False)

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
