#!/usr/bin/env python3
"""Train + test + interpret one LightGBM resistance classifier per antibiotic.

Importance is mean |SHAP| over the training rows (TreeExplainer: exact and
cheap for trees), reported alongside the standard SHAP beeswarm plot.
"""
import numpy as np
import shap
from lightgbm import LGBMClassifier

import common
from common import plt

MODEL_NAME = "lightgbm"

DEFAULTS = {
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
    # minority class contributes equally to training. More stable than manual
    # scale_pos_weight on tiny, skewed AMR splits.
    "is_unbalance": True,
    "random_state": 42,
    "verbose": -1,
}


def main():
    args = common.train_parser("Train a LightGBM resistance model").parse_args()
    common.prepare_outputs(args)

    genes = common.read_genes(args.train_features, args.all_antibiotics)
    x_train, y_train, _ = common.load_xy(args.train_features, args.antibiotic, genes)
    x_test, y_test, test_runs = common.load_xy(args.test_features, args.antibiotic, genes)

    base = {"model": MODEL_NAME, "antibiotic": args.antibiotic,
            "n_train": int(len(x_train)), "n_test": int(len(x_test)),
            "n_features": len(genes)}

    if not common.enough_to_fit(y_train, len(x_test)):
        common.write_skipped(args, base, MODEL_NAME)
        return

    n_pos = int(y_train.sum())
    n_neg = int(len(y_train) - n_pos)

    hyperparams = {**DEFAULTS, **common.load_tuned(args.params_input)}

    model = LGBMClassifier(**hyperparams)
    model.fit(x_train, y_train)
    model.booster_.save_model(args.model_output)
    common.write_json(args.params_output,
                      {**base, "hyperparameters": hyperparams,
                       "train_class_balance": {"R": n_pos, "S": n_neg}})

    y_proba = model.predict_proba(x_test)[:, 1]
    y_pred = model.predict(x_test)
    metrics = common.compute_metrics(base, y_test, y_pred, y_proba)
    common.write_json(args.metrics_output, metrics)
    common.write_predictions(args.predictions_output, test_runs, y_test, y_pred, y_proba)

    shap_values = shap.TreeExplainer(model).shap_values(x_train)
    # LightGBM binary TreeExplainer may return a list [class0, class1]; take positive.
    if isinstance(shap_values, list):
        shap_values = shap_values[1]
    shap_values = np.asarray(shap_values)
    common.write_importance(args.importance_output, genes,
                            np.abs(shap_values).mean(axis=0))

    shap.summary_plot(shap_values, x_train, show=False, max_display=20)
    plt.title(f"{args.antibiotic} ({MODEL_NAME}) SHAP")
    plt.savefig(args.importance_plot_output, dpi=150, bbox_inches="tight")
    plt.close()

    common.print_json(metrics)


if __name__ == "__main__":
    main()
