#!/usr/bin/env python3
# /// script
# dependencies = ["pandas", "numpy", "scikit-learn", "matplotlib", "torch"]
#
# [[tool.uv.index]]
# name = "pytorch-cpu"
# url = "https://download.pytorch.org/whl/cpu"
# explicit = true
#
# [tool.uv.sources]
# torch = { index = "pytorch-cpu" }
# ///
"""Train + test + interpret an MLP resistance classifier per antibiotic.

A fully-connected network on the AMR-gene presence/absence matrix -- the
principled neural architecture for tabular features (no spatial structure
to convolve, so a CNN would be unmotivated). Same features as the tree
models. Produces: saved model (torch), hyperparameters, test predictions,
metrics, and permutation feature importance (CSV + PNG). CPU only.
"""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
)

MODEL_NAME = "nn"
torch.manual_seed(42)

HYPERPARAMS = {"epochs": 200, "batch_size": 32, "learning_rate": 1e-3,
               "hidden": [128, 64], "dropout": 0.3, "weight_decay": 1e-4}


class MLP(nn.Module):
    def __init__(self, input_size, hidden, dropout):
        super().__init__()
        layers, prev = [], input_size
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(dropout)]
            prev = h
        layers += [nn.Linear(prev, 2)]
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


def gene_columns(df, all_antibiotics):
    non_gene = {"run", "collection_date", "year", "location", *all_antibiotics}
    return [c for c in df.columns if c not in non_gene]


def load_xy(csv_path, antibiotic, genes):
    df = pd.read_csv(csv_path).dropna(subset=[antibiotic])
    y = (df[antibiotic] == "R").astype(int).values
    x = df[genes].fillna(0).values.astype("float32")
    return x, y, df["run"]


def train(model, x, y, hp, class_weights):
    opt = torch.optim.Adam(model.parameters(), lr=hp["learning_rate"],
                           weight_decay=hp["weight_decay"])
    criterion = nn.CrossEntropyLoss(weight=class_weights)
    xt, yt = torch.tensor(x), torch.tensor(y)
    n, bs = len(x), min(hp["batch_size"], len(x))
    model.train()
    for _ in range(hp["epochs"]):
        perm = torch.randperm(n)
        for i in range(0, n, bs):
            idx = perm[i:i + bs]
            opt.zero_grad()
            loss = criterion(model(xt[idx]), yt[idx])
            loss.backward()
            opt.step()


@torch.no_grad()
def predict_proba(model, x):
    model.eval()
    return torch.softmax(model(torch.tensor(x)), dim=1)[:, 1].numpy()


def placeholder_png(path, text):
    fig, ax = plt.subplots(figsize=(6, 2))
    ax.text(0.5, 0.5, text, ha="center", va="center", wrap=True)
    ax.axis("off")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    p = argparse.ArgumentParser(description="Train an MLP resistance model")
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
    return p.parse_args()


def main():
    args = parse_args()
    for out in (args.model_output, args.params_output, args.metrics_output,
                args.predictions_output, args.importance_output, args.importance_plot_output):
        Path(out).parent.mkdir(parents=True, exist_ok=True)

    genes = gene_columns(pd.read_csv(args.train_features), args.all_antibiotics)
    x_train, y_train, _ = load_xy(args.train_features, args.antibiotic, genes)
    x_test, y_test, test_runs = load_xy(args.test_features, args.antibiotic, genes)

    base = {"model": MODEL_NAME, "antibiotic": args.antibiotic,
            "n_train": int(len(x_train)), "n_test": int(len(x_test)),
            "n_features": len(genes)}

    if len(x_train) < 2 or len(x_test) == 0 or len(np.unique(y_train)) < 2:
        base["skipped"] = ("need both R and S (>=2) in train and >=1 test sample "
                           "(not enough labeled data for this antibiotic)")
        Path(args.model_output).write_bytes(b"")
        json.dump(base, open(args.params_output, "w"), indent=2)
        json.dump(base, open(args.metrics_output, "w"), indent=2)
        pd.DataFrame(columns=["run", "actual", "predicted", "probability_resistant"]
                     ).to_csv(args.predictions_output, index=False)
        pd.DataFrame(columns=["gene", "importance"]).to_csv(args.importance_output, index=False)
        placeholder_png(args.importance_plot_output, f"{args.antibiotic} ({MODEL_NAME}): skipped\n{base['skipped']}")
        print(json.dumps(base, indent=2))
        return

    n_pos = int(y_train.sum())
    n_neg = int(len(y_train) - n_pos)
    class_weights = torch.tensor([len(y_train) / (2 * max(n_neg, 1)),
                                  len(y_train) / (2 * max(n_pos, 1))], dtype=torch.float32)

    model = MLP(len(genes), HYPERPARAMS["hidden"], HYPERPARAMS["dropout"])
    train(model, x_train, y_train, HYPERPARAMS, class_weights)

    torch.save({"state_dict": model.state_dict(), "input_size": len(genes),
                "genes": genes, "hyperparameters": HYPERPARAMS}, args.model_output)
    json.dump({**base, "hyperparameters": HYPERPARAMS,
               "train_class_balance": {"R": n_pos, "S": n_neg}},
              open(args.params_output, "w"), indent=2)

    y_proba = predict_proba(model, x_test)
    y_pred = (y_proba >= 0.5).astype(int)
    metrics = {**base,
               "accuracy": accuracy_score(y_test, y_pred),
               "balanced_accuracy": balanced_accuracy_score(y_test, y_pred),
               "f1": f1_score(y_test, y_pred, zero_division=0),
               "roc_auc": roc_auc_score(y_test, y_proba) if len(np.unique(y_test)) > 1 else None}
    json.dump(metrics, open(args.metrics_output, "w"), indent=2)
    pd.DataFrame({"run": test_runs, "actual": y_test,
                  "predicted": y_pred, "probability_resistant": y_proba}
                 ).to_csv(args.predictions_output, index=False)

    # Permutation importance (model-agnostic): drop in held-out score when
    # each gene column is shuffled.
    rng = np.random.default_rng(42)
    scored = lambda pr: (roc_auc_score(y_test, pr) if len(np.unique(y_test)) > 1
                         else accuracy_score(y_test, (pr >= 0.5).astype(int)))
    base_score = scored(y_proba)
    drops = []
    for j in range(len(genes)):
        xp = x_test.copy()
        rng.shuffle(xp[:, j])
        drops.append(base_score - scored(predict_proba(model, xp)))
    importance = (pd.DataFrame({"gene": genes, "importance": drops})
                  .sort_values("importance", ascending=False).reset_index(drop=True))
    importance.to_csv(args.importance_output, index=False)

    top = importance.head(20)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(top["gene"][::-1], top["importance"][::-1])
    ax.set_xlabel("permutation importance (drop in test score)")
    ax.set_title(f"{args.antibiotic} ({MODEL_NAME}) feature importance")
    fig.savefig(args.importance_plot_output, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
