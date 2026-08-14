#!/usr/bin/env python3
"""Optuna hyperparameter search for LightGBM resistance classifiers.

Searches the training split with stratified k-fold CV and writes the winning
hyperparameters for train_lightgbm.py to pick up via --params-input. The driver
lives in common.run_tuning; this file is the search space and the fit step.
"""
from lightgbm import LGBMClassifier

import common


def suggest(trial):
    return {
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
        # Matches the trainer's imbalance handling so CV scores reflect how the
        # model is actually fitted.
        "is_unbalance": True,
        "random_state": 42,
        "verbose": -1,
    }


def fit_predict(hyperparams, x_train, y_train, x_val):
    model = LGBMClassifier(**hyperparams)
    model.fit(x_train, y_train)
    return model.predict_proba(x_val)[:, 1]


if __name__ == "__main__":
    args = common.tune_parser("Tune LightGBM hyperparameters with Optuna", 30).parse_args()
    common.run_tuning(args, suggest, fit_predict)
