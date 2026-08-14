#!/usr/bin/env python3
"""Optuna hyperparameter search for the MLP resistance classifier.

Searches the training split with stratified k-fold CV and writes the winning
hyperparameters for train_nn.py to pick up via --params-input. The driver lives
in common.run_tuning; this file is the search space and the fit step.

Note `hidden` is assembled from two scalar suggestions into the layer-size list
train_nn.py expects. common.run_tuning replays this same function over the best
trial, so the list -- not the loose hidden0/hidden1 scalars -- is what gets
written out.
"""
import numpy as np
import torch

import common
import mlp

torch.manual_seed(42)
np.random.seed(42)


def suggest(trial):
    return {
        "epochs": trial.suggest_int("epochs", 50, 300, log=True),
        "batch_size": trial.suggest_categorical("batch_size", [8, 16, 32, 64]),
        "learning_rate": trial.suggest_float("learning_rate", 1e-4, 1e-2, log=True),
        "weight_decay": trial.suggest_float("weight_decay", 1e-6, 1e-3, log=True),
        "dropout": trial.suggest_float("dropout", 0.1, 0.6),
        "hidden": [trial.suggest_int("hidden0", 32, 256, log=True),
                   trial.suggest_int("hidden1", 16, 128, log=True)],
    }


def fit_predict(hyperparams, x_train, y_train, x_val):
    return mlp.fit_predict(hyperparams, x_train, y_train, x_val,
                           common.class_weights(y_train))


if __name__ == "__main__":
    args = common.tune_parser("Tune MLP hyperparameters with Optuna", 20).parse_args()
    common.run_tuning(args, suggest, fit_predict)
