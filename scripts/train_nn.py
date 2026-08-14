#!/usr/bin/env python3
"""Train + test + interpret an MLP resistance classifier per antibiotic.

Same features as the tree models; the architecture and training loop live in
mlp.py, shared with tune_nn.py. Importance is permutation importance on the
test split -- model-agnostic, and the natural choice for a net where SHAP's
tree fast path does not apply.
"""
import numpy as np
import torch
from sklearn.metrics import accuracy_score, roc_auc_score

import common
import mlp

MODEL_NAME = "nn"
torch.manual_seed(42)


def permutation_importance(model, x_test, y_test, baseline_proba, rng):
    """Drop in held-out score when each gene column is shuffled."""
    def score(proba):
        if len(np.unique(y_test)) > 1:
            return roc_auc_score(y_test, proba)
        return accuracy_score(y_test, (proba >= 0.5).astype(int))

    base_score = score(baseline_proba)
    drops = []
    for j in range(x_test.shape[1]):
        permuted = x_test.copy()
        rng.shuffle(permuted[:, j])
        drops.append(base_score - score(mlp.predict_proba(model, permuted)))
    return drops


def main():
    args = common.train_parser("Train an MLP resistance model").parse_args()
    common.prepare_outputs(args)

    genes = common.read_genes(args.train_features, args.all_antibiotics)
    x_train, y_train, _ = common.load_arrays(args.train_features, args.antibiotic, genes)
    x_test, y_test, test_runs = common.load_arrays(args.test_features, args.antibiotic, genes)

    base = {"model": MODEL_NAME, "antibiotic": args.antibiotic,
            "n_train": int(len(x_train)), "n_test": int(len(x_test)),
            "n_features": len(genes)}

    # Two training rows minimum: a single row cannot carry both classes.
    if not common.enough_to_fit(y_train, len(x_test), min_train=2):
        common.write_skipped(args, base, MODEL_NAME)
        return

    n_pos = int(y_train.sum())
    n_neg = int(len(y_train) - n_pos)
    hyperparams = {**mlp.DEFAULTS, **common.load_tuned(args.params_input)}

    model = mlp.MLP(len(genes), hyperparams["hidden"], hyperparams["dropout"])
    mlp.train(model, x_train, y_train, hyperparams, common.class_weights(y_train))

    torch.save({"state_dict": model.state_dict(), "input_size": len(genes),
                "genes": genes, "hyperparameters": hyperparams}, args.model_output)
    common.write_json(args.params_output,
                      {**base, "hyperparameters": hyperparams,
                       "train_class_balance": {"R": n_pos, "S": n_neg}})

    y_proba = mlp.predict_proba(model, x_test)
    y_pred = (y_proba >= 0.5).astype(int)
    metrics = common.compute_metrics(base, y_test, y_pred, y_proba)
    common.write_json(args.metrics_output, metrics)
    common.write_predictions(args.predictions_output, test_runs, y_test, y_pred, y_proba)

    drops = permutation_importance(model, x_test, y_test, y_proba,
                                   np.random.default_rng(42))
    importance = common.write_importance(args.importance_output, genes, drops)
    common.importance_plot(args.importance_plot_output, importance,
                           f"{args.antibiotic} ({MODEL_NAME}) feature importance",
                           "permutation importance (drop in test score)")

    common.print_json(metrics)


if __name__ == "__main__":
    main()
