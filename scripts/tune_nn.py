#!/usr/bin/env python3
"""Optuna hyperparameter search for the MLP resistance classifier.

Searches on the training split using stratified k-fold CV and writes the best
hyperparameters as JSON. The training script loads this JSON via --params-input.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import optuna
import pandas as pd
import torch
import torch.nn as nn
from sklearn.model_selection import StratifiedKFold

torch.manual_seed(42)
np.random.seed(42)
optuna.logging.set_verbosity(optuna.logging.WARNING)


def gene_columns(df, all_antibiotics):
    non_gene = {"run", "collection_date", "year", "location", *all_antibiotics}
    return [c for c in df.columns if c not in non_gene]


def load_xy(csv_path, antibiotic, genes):
    df = pd.read_csv(csv_path).dropna(subset=[antibiotic])
    y = (df[antibiotic] == "R").astype(int).values
    x = df[genes].fillna(0).values.astype("float32")
    return x, y


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


def train_mlp(model, x, y, hp, class_weights):
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


def class_weights_tensor(y):
    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    return torch.tensor([len(y) / (2 * max(n_neg, 1)),
                         len(y) / (2 * max(n_pos, 1))], dtype=torch.float32)


@torch.no_grad()
def predict_proba(model, x):
    model.eval()
    return torch.softmax(model(torch.tensor(x)), dim=1)[:, 1].numpy()


def cv_score(x, y, hyperparams, n_splits=3, seed=42):
    if len(np.unique(y)) < 2:
        return 0.0
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    aucs = []
    for tr_idx, val_idx in skf.split(x, y):
        x_tr, x_val = x[tr_idx], x[val_idx]
        y_tr, y_val = y[tr_idx], y[val_idx]
        model = MLP(x.shape[1], hyperparams["hidden"], hyperparams["dropout"])
        train_mlp(model, x_tr, y_tr, hyperparams, class_weights_tensor(y_tr))
        proba = predict_proba(model, x_val)
        if len(np.unique(y_val)) > 1:
            from sklearn.metrics import roc_auc_score
            aucs.append(roc_auc_score(y_val, proba))
    return float(np.mean(aucs)) if aucs else 0.0


def objective(trial, x, y):
    hyperparams = {
        "epochs": trial.suggest_int("epochs", 50, 300, log=True),
        "batch_size": trial.suggest_categorical("batch_size", [8, 16, 32, 64]),
        "learning_rate": trial.suggest_float("learning_rate", 1e-4, 1e-2, log=True),
        "weight_decay": trial.suggest_float("weight_decay", 1e-6, 1e-3, log=True),
        "dropout": trial.suggest_float("dropout", 0.1, 0.6),
        "hidden": [trial.suggest_int("hidden0", 32, 256, log=True),
                   trial.suggest_int("hidden1", 16, 128, log=True)],
    }
    return cv_score(x, y, hyperparams)


def parse_args():
    p = argparse.ArgumentParser(description="Tune MLP hyperparameters with Optuna")
    p.add_argument("--train-features", required=True)
    p.add_argument("--antibiotic", required=True)
    p.add_argument("--all-antibiotics", nargs="+", required=True)
    p.add_argument("--output", required=True, help="JSON file for best hyperparameters")
    p.add_argument("--n-trials", type=int, default=20)
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
