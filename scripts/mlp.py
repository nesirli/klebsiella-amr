#!/usr/bin/env python3
"""The MLP shared by train_nn.py and tune_nn.py.

A fully-connected net over the AMR-gene presence/absence matrix -- the
principled neural architecture for tabular features (no spatial structure to
convolve, so a CNN would be unmotivated). CPU only.
"""
import numpy as np
import torch
import torch.nn as nn

DEFAULTS = {"epochs": 200, "batch_size": 32, "learning_rate": 1e-3,
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


def train(model, x, y, hyperparams, class_weights):
    opt = torch.optim.Adam(model.parameters(), lr=hyperparams["learning_rate"],
                           weight_decay=hyperparams["weight_decay"])
    criterion = nn.CrossEntropyLoss(weight=class_weights)
    xt, yt = torch.tensor(x), torch.tensor(y)
    n, bs = len(x), min(hyperparams["batch_size"], len(x))
    model.train()
    for _ in range(hyperparams["epochs"]):
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


def fit_predict(hyperparams, x_train, y_train, x_val, class_weights):
    """Fresh model, trained on (x_train, y_train) -> P(resistant) for x_val."""
    model = MLP(x_train.shape[1], hyperparams["hidden"], hyperparams["dropout"])
    train(model, x_train, y_train, hyperparams, class_weights)
    return predict_proba(model, np.asarray(x_val, dtype="float32"))
