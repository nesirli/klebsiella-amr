#!/usr/bin/env python3
"""Optuna search for DNABERT-2 -- the sequence model's counterpart to the
tabular tuners.

This cannot use common.run_tuning: that driver loads the feature matrix as a
float array, while this model consumes the nucleotide strings themselves. It
reuses the same fold logic (common.cv_roc_auc) and writes the same
{"hyperparameters", "tuning"} handoff, so the tuner -> trainer contract is
identical to the other three models.

Fine-tuning a 117M-parameter encoder costs orders of magnitude more per fit
than fitting a GBM, so the search is deliberately small and runs on a capped,
class-balanced subsample of the training split. The goal is a sane learning
rate and epoch count, not an exhaustive sweep: a fold here is minutes, not
milliseconds. Anything the search does not name keeps the trainer's default.
"""
import numpy as np
import optuna
import torch
from transformers import AutoModel, AutoTokenizer

import common
import train_dnabert as td

# Every fold refits from the pretrained checkpoint, so a trial costs
# n_splits full fine-tunes. Keep both small.
DEFAULT_TRIALS = 5
DEFAULT_MAX_TRAIN = 120


def suggest(trial):
    """One complete hyperparameter set. encode_batch is a memory knob, not a
    quality one, so it stays fixed and out of the search."""
    return {
        "learning_rate": trial.suggest_float("learning_rate", 5e-6, 5e-5, log=True),
        "epochs": trial.suggest_int("epochs", 1, 3),
        "batch_size": trial.suggest_categorical("batch_size", [2, 4]),
        "max_length": trial.suggest_categorical("max_length", [64, 128]),
        "weight_decay": trial.suggest_float("weight_decay", 1e-4, 1e-1, log=True),
        "encode_batch": 8,
    }


def fit_predict(hyperparams, x_train, y_train, x_val, tokenizer):
    """Fine-tune from the pretrained checkpoint and score the held-out fold.

    The encoder is reloaded per fold on purpose: carrying weights across folds
    would leak the validation half of one fold into the next.
    """
    encoder = AutoModel.from_pretrained(td.BASE_MODEL, trust_remote_code=True)
    td.force_pytorch_attention()
    model = td.DNABertGeneClassifier(encoder, encoder.config.hidden_size).to(td.DEVICE)

    td.train(model, tokenizer, list(x_train), y_train, hyperparams,
             common.class_weights(y_train))

    model.eval()
    embeddings, _ = td.sample_embeddings(model, tokenizer, list(x_val), hyperparams)
    pooled = torch.stack([e.mean(0) for e in embeddings]).to(td.DEVICE)
    return td.proba_from_pooled(model, pooled)


def balanced_subsample(x, y, cap, seed=42):
    """At most `cap` samples, keeping both classes proportionally represented."""
    if cap <= 0 or cap >= len(y):
        return x, y
    rng = np.random.default_rng(seed)
    keep = []
    for label in np.unique(y):
        idx = np.flatnonzero(y == label)
        take = max(1, round(cap * len(idx) / len(y)))
        keep.append(rng.choice(idx, size=min(take, len(idx)), replace=False))
    picked = np.sort(np.concatenate(keep))
    return x[picked], y[picked]


def parse_args():
    p = common.tune_parser("Tune DNABERT-2 hyperparameters", DEFAULT_TRIALS)
    p.add_argument("--max-train", type=int, default=DEFAULT_MAX_TRAIN,
                   help="cap the tuning subsample; -1 uses the whole split")
    return p.parse_args()


def main():
    args = parse_args()
    optuna.logging.set_verbosity(optuna.logging.WARNING)

    genes = common.read_genes(args.train_features, args.all_antibiotics)
    train_df = common.load_labelled(args.train_features, args.antibiotic,
                                    dtype={g: str for g in genes})
    y = (train_df[args.antibiotic] == "R").astype(int).to_numpy()

    # Object array, not a numeric one: each element is a sample's list of
    # (gene, sequence) pairs. Built empty-then-filled so numpy keeps it 1-D
    # instead of trying to broadcast the ragged lists into a matrix.
    gene_lists = td.sample_gene_lists(train_df, genes)
    x = np.empty(len(gene_lists), dtype=object)
    x[:] = gene_lists

    if (len(np.unique(y)) < 2 or len(y) < args.n_splits * 2
            or np.bincount(y).min() < 2):
        common.write_tuning_result(args.output, {}, {"skipped": common.TUNE_SKIP_REASON})
        return

    x, y = balanced_subsample(x, y, args.max_train, args.seed)

    tokenizer = AutoTokenizer.from_pretrained(td.BASE_MODEL, trust_remote_code=True)

    def objective(trial):
        hyperparams = suggest(trial)
        return common.cv_roc_auc(
            x, y,
            lambda xt, yt, xv: fit_predict(hyperparams, xt, yt, xv, tokenizer),
            args.n_splits, args.seed,
        )

    study = optuna.create_study(direction="maximize",
                                sampler=optuna.samplers.TPESampler(seed=args.seed))
    study.optimize(objective, n_trials=args.n_trials)

    common.write_tuning_result(
        args.output,
        suggest(optuna.trial.FixedTrial(study.best_params)),
        {"n_trials": args.n_trials, "n_splits": args.n_splits,
         "best_cv_roc_auc": float(study.best_value),
         "tuned_on_samples": int(len(y)), "device": td.DEVICE.type},
    )


if __name__ == "__main__":
    main()
