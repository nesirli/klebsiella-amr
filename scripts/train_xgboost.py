#!/usr/bin/env python3
"""Train + test + interpret one XGBoost resistance classifier per antibiotic.

Importance is mean |SHAP| over the training rows (TreeExplainer: exact and
cheap for trees), reported alongside the standard SHAP beeswarm plot.
"""
import numpy as np
import shap
from xgboost import XGBClassifier

import common
from common import plt

MODEL_NAME = "xgboost"

DEFAULTS = {
    "objective": "binary:logistic",
    "eval_metric": ["logloss", "auc"],
    "n_estimators": 200,
    "max_depth": 4,
    "learning_rate": 0.1,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "min_child_weight": 1,
    "max_delta_step": 1,
    "reg_alpha": 0.1,
    "reg_lambda": 0.1,
    "base_score": 0.5,
    "random_state": 42,
}


def main():
    args = common.train_parser("Train an XGBoost resistance model").parse_args()
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

    # Class imbalance: scale_pos_weight is data-dependent, so it is computed
    # here rather than searched, and a tuned value would override it.
    hyperparams = {**DEFAULTS, "scale_pos_weight": (n_neg / n_pos) if n_pos else 1.0}
    hyperparams.update(common.load_tuned(args.params_input))

    model = XGBClassifier(**hyperparams)
    model.fit(x_train, y_train)
    model.get_booster().save_model(args.model_output)
    common.write_json(args.params_output,
                      {**base, "hyperparameters": hyperparams,
                       "train_class_balance": {"R": n_pos, "S": n_neg}})

    y_proba = model.predict_proba(x_test)[:, 1]
    y_pred = model.predict(x_test)
    metrics = common.compute_metrics(base, y_test, y_pred, y_proba)
    common.write_json(args.metrics_output, metrics)
    common.write_predictions(args.predictions_output, test_runs, y_test, y_pred, y_proba)

    shap_values = np.asarray(shap.TreeExplainer(model).shap_values(x_train))
    common.write_importance(args.importance_output, genes,
                            np.abs(shap_values).mean(axis=0))

    shap.summary_plot(shap_values, x_train, show=False, max_display=20)
    plt.title(f"{args.antibiotic} ({MODEL_NAME}) SHAP")
    plt.savefig(args.importance_plot_output, dpi=150, bbox_inches="tight")
    plt.close()

    common.print_json(metrics)


if __name__ == "__main__":
    main()
