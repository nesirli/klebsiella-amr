#!/usr/bin/env python3
"""Optuna hyperparameter search for XGBoost resistance classifiers.

Searches the training split with stratified k-fold CV and writes the winning
hyperparameters for train_xgboost.py to pick up via --params-input. The driver
lives in common.run_tuning; this file is the search space and the fit step.
"""
from xgboost import XGBClassifier

import common


def suggest(trial):
    return {
        "objective": "binary:logistic",
        "eval_metric": "logloss",
        "n_estimators": trial.suggest_int("n_estimators", 100, 500, log=True),
        "max_depth": trial.suggest_int("max_depth", 2, 8),
        "learning_rate": trial.suggest_float("learning_rate", 1e-3, 0.3, log=True),
        "subsample": trial.suggest_float("subsample", 0.5, 1.0),
        "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, 1.0),
        "min_child_weight": trial.suggest_int("min_child_weight", 1, 10),
        "max_delta_step": trial.suggest_int("max_delta_step", 0, 10),
        "reg_alpha": trial.suggest_float("reg_alpha", 1e-8, 10.0, log=True),
        "reg_lambda": trial.suggest_float("reg_lambda", 1e-8, 10.0, log=True),
        "random_state": 42,
    }


def fit_predict(hyperparams, x_train, y_train, x_val):
    # scale_pos_weight is a property of the fold, not of the search space, so
    # it is recomputed per fold and left out of the written hyperparameters.
    n_pos = int(y_train.sum())
    n_neg = int(len(y_train) - n_pos)
    model = XGBClassifier(**hyperparams,
                          scale_pos_weight=(n_neg / n_pos) if n_pos else 1.0)
    model.fit(x_train, y_train)
    return model.predict_proba(x_val)[:, 1]


if __name__ == "__main__":
    args = common.tune_parser("Tune XGBoost hyperparameters with Optuna", 30).parse_args()
    common.run_tuning(args, suggest, fit_predict)
