import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ANTIBIOTICS = ["amikacin", "ciprofloxacin", "ceftazidime", "meropenem"]

# script, model file extension, number of tuning trials
MODELS = {
    "xgboost": ("scripts/train_xgboost.py", "scripts/tune_xgboost.py", "json", 3),
    "lightgbm": ("scripts/train_lightgbm.py", "scripts/tune_lightgbm.py", "txt", 3),
    "nn": ("scripts/train_nn.py", "scripts/tune_nn.py", "pt", 2),
}


def _make_feature_csv(path, n_samples=50, n_genes=20, seed=42):
    rng = np.random.default_rng(seed)
    genes = [f"gene_{i}" for i in range(n_genes)]
    rows = []
    for i in range(n_samples):
        row = {
            "run": f"SRR{i:05d}",
            "collection_date": "2020-01-01",
            "year": 2020,
            "location": "Test",
        }
        for abx in ANTIBIOTICS:
            row[abx] = "R" if rng.random() > 0.5 else "S"
        for g in genes:
            row[g] = int(rng.random() > 0.5)
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)


def _run(*args):
    subprocess.run([sys.executable, *[str(a) for a in args]], check=True)


def _train(model, abx, train, test, out_dir, params_input=None):
    """Run one training script. Every output goes under out_dir.

    Output paths are built here rather than by the caller so that no artifact
    can escape the temporary directory and land in the repository root.
    """
    script, _, ext, _ = MODELS[model]
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_name = f"{abx}_shap.png" if model in ("xgboost", "lightgbm") else f"{abx}_importance.png"
    args = [
        script,
        "--train-features", train,
        "--test-features", test,
        "--antibiotic", abx,
        "--all-antibiotics", *ANTIBIOTICS,
        "--model-output", out_dir / f"{abx}_model.{ext}",
        "--params-output", out_dir / f"{abx}_params.json",
        "--metrics-output", out_dir / f"{abx}_metrics.json",
        "--predictions-output", out_dir / f"{abx}_predictions.csv",
        "--importance-output", out_dir / f"{abx}_importance.csv",
        "--importance-plot-output", out_dir / plot_name,
    ]
    if params_input:
        args += ["--params-input", params_input]
    _run(*args)
    return json.loads((out_dir / f"{abx}_metrics.json").read_text())


def _tune(model, abx, train, out_path, n_splits=3):
    _, script, _, n_trials = MODELS[model]
    _run(script,
         "--train-features", train,
         "--antibiotic", abx,
         "--all-antibiotics", *ANTIBIOTICS,
         "--output", out_path,
         "--n-trials", n_trials,
         "--n-splits", n_splits)
    return json.loads(Path(out_path).read_text())


@pytest.fixture
def features(tmp_path):
    train, test = tmp_path / "train.csv", tmp_path / "test.csv"
    _make_feature_csv(train, n_samples=80)
    _make_feature_csv(test, n_samples=20, seed=1)
    return train, test


@pytest.mark.parametrize("model", list(MODELS))
def test_model_runs(model, features, tmp_path):
    train, test = features
    for abx in ANTIBIOTICS:
        metrics = _train(model, abx, train, test, tmp_path / model)
        for key in ("f1", "precision", "recall", "pr_auc", "confusion_matrix"):
            assert key in metrics
        assert metrics["n_train"] == 80
        assert metrics["n_test"] == 20


@pytest.mark.parametrize("model", list(MODELS))
def test_model_emits_every_artifact(model, features, tmp_path):
    """make treats a missing output as a failed target, so all six must exist."""
    train, test = features
    out_dir = tmp_path / model
    _train(model, "amikacin", train, test, out_dir)
    ext = MODELS[model][2]
    plot = "amikacin_shap.png" if model in ("xgboost", "lightgbm") else "amikacin_importance.png"
    for name in (f"amikacin_model.{ext}", "amikacin_params.json", "amikacin_metrics.json",
                 "amikacin_predictions.csv", "amikacin_importance.csv", plot):
        assert (out_dir / name).is_file(), f"{model} did not write {name}"

    importance = pd.read_csv(out_dir / "amikacin_importance.csv")
    assert list(importance.columns) == ["gene", "importance"]
    assert importance["importance"].is_monotonic_decreasing


@pytest.mark.parametrize("model", list(MODELS))
def test_tuning_pipeline(model, features, tmp_path):
    """End-to-end: tune hyperparameters, then train with the tuned params."""
    train, test = features
    out_dir = tmp_path / model
    out_dir.mkdir()

    for abx in ANTIBIOTICS:
        tune_out = out_dir / f"tune_{abx}.json"
        result = _tune(model, abx, train, tune_out)
        assert "best_cv_roc_auc" in result["tuning"] or "skipped" in result["tuning"]

        metrics = _train(model, abx, train, test, out_dir, params_input=tune_out)
        assert "f1" in metrics


def test_tuning_result_separates_hyperparameters_from_bookkeeping(features, tmp_path):
    """Search bookkeeping must not reach the estimator constructor.

    n_trials/n_splits/best_cv_roc_auc once sat beside the real hyperparameters
    in the same flat dict, and the trainer splatted the lot into the model.
    """
    train, _ = features
    result = _tune("xgboost", "amikacin", train, tmp_path / "tune.json")

    assert set(result) == {"hyperparameters", "tuning"}
    assert set(result["tuning"]) == {"n_trials", "n_splits", "best_cv_roc_auc"}
    assert not set(result["hyperparameters"]) & set(result["tuning"])


def test_tuned_nn_architecture_reaches_the_model(features, tmp_path):
    """The tuned layer sizes must survive the tuner -> trainer handoff.

    Optuna suggests hidden0/hidden1 as scalars; the trainer wants a `hidden`
    list. When the handoff passed the flat scalars through, the trainer
    silently kept its default architecture and the search was wasted.
    """
    train, test = features
    tuned = _tune("nn", "amikacin", train, tmp_path / "tune.json", n_splits=2)

    hidden = tuned["hyperparameters"]["hidden"]
    assert isinstance(hidden, list) and len(hidden) == 2
    assert "hidden0" not in tuned["hyperparameters"]

    out_dir = tmp_path / "nn"
    _train("nn", "amikacin", train, test, out_dir, params_input=tmp_path / "tune.json")

    params = json.loads((out_dir / "amikacin_params.json").read_text())
    assert params["hyperparameters"]["hidden"] == hidden

    import torch
    saved = torch.load(out_dir / "amikacin_model.pt", weights_only=False)
    layer_widths = [t.shape[0] for name, t in saved["state_dict"].items()
                    if name.endswith("weight")]
    assert layer_widths[:2] == hidden, "saved net does not use the tuned widths"


def test_skipped_tuning_does_not_poison_training(features, tmp_path):
    """A tuner that gave up must not hand the trainer junk to splat."""
    train, test = features
    single_class = pd.read_csv(train)
    single_class["amikacin"] = "R"
    flat = tmp_path / "flat.csv"
    single_class.to_csv(flat, index=False)

    tuned = _tune("xgboost", "amikacin", flat, tmp_path / "tune.json")
    assert "skipped" in tuned["tuning"]
    assert tuned["hyperparameters"] == {}

    # Training on the real (two-class) split still works with that file.
    metrics = _train("xgboost", "amikacin", train, test, tmp_path / "xgb",
                     params_input=tmp_path / "tune.json")
    assert "f1" in metrics


@pytest.mark.parametrize("model", list(MODELS))
def test_single_class_training_data_is_skipped_cleanly(model, features, tmp_path):
    """Too little signal to fit: still emit every artifact, marked skipped."""
    train, test = features
    single_class = pd.read_csv(train)
    single_class["amikacin"] = "R"
    flat = tmp_path / "flat.csv"
    single_class.to_csv(flat, index=False)

    out_dir = tmp_path / model
    metrics = _train(model, "amikacin", flat, test, out_dir)
    assert "skipped" in metrics
    assert (out_dir / "amikacin_predictions.csv").is_file()
    assert (out_dir / "amikacin_importance.csv").is_file()
