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


def test_threshold_rescues_a_ranking_the_default_cut_would_waste():
    """A perfect ranking under 0.5 must still produce positive calls.

    DNABERT-2 on amikacin ranked the test set at ROC AUC 0.980 and then
    labelled all 178 isolates susceptible: its probabilities never reached the
    hardcoded 0.5, because training was 74 resistant against 228 susceptible.
    """
    import sys
    sys.path.insert(0, "scripts")
    import common

    y = np.array([0] * 20 + [1] * 5)
    # Perfectly separable, but everything sits well below 0.5.
    proba = np.concatenate([np.linspace(0.01, 0.05, 20), np.linspace(0.10, 0.14, 5)])

    assert (proba >= 0.5).sum() == 0, "fixture should defeat the fixed cut"
    threshold = common.choose_threshold(y, proba)
    predicted = (proba >= threshold).astype(int)
    assert predicted.sum() == 5 and (predicted == y).all()


def test_threshold_falls_back_when_one_class_is_absent():
    import sys
    sys.path.insert(0, "scripts")
    import common

    assert common.choose_threshold(np.zeros(10), np.linspace(0, 1, 10)) == 0.5


def test_dnabert_tuner_writes_the_same_handoff_as_the_others(tmp_path):
    """DNABERT-2 has its own driver (it consumes sequences, not a float
    matrix), so its handoff file has to be checked against the shared contract
    independently. Uses the single-class skip path: exercising a real trial
    would fine-tune a 117M-parameter encoder.
    """
    rows = []
    for i in range(6):
        row = {"run": f"SRR{i:05d}", "collection_date": "2020-01-01",
               "year": 2020, "location": "Test"}
        for abx in ANTIBIOTICS:
            row[abx] = "R"  # single class -> tuner must skip, not crash
        row["blaKPC-3"] = "ATGCGT" * 10
        rows.append(row)
    seqs = tmp_path / "train_sequences.csv"
    pd.DataFrame(rows).to_csv(seqs, index=False)

    out = tmp_path / "tune_amikacin.json"
    _run("scripts/tune_dnabert.py",
         "--train-features", seqs,
         "--antibiotic", "amikacin",
         "--all-antibiotics", *ANTIBIOTICS,
         "--output", out,
         "--n-trials", 1,
         "--n-splits", 2)

    result = json.loads(out.read_text())
    assert set(result) == {"hyperparameters", "tuning"}
    assert "skipped" in result["tuning"]
    assert result["hyperparameters"] == {}


def test_dnabert_search_space_covers_the_trainer_knobs():
    """A tuned key the trainer ignores is a silently wasted search."""
    import sys
    sys.path.insert(0, "scripts")
    import optuna

    import tune_dnabert

    params = tune_dnabert.suggest(optuna.trial.FixedTrial({
        "learning_rate": 2e-5, "epochs": 2, "batch_size": 2,
        "max_length": 128, "weight_decay": 1e-2,
    }))
    # Every searched key must be one train_dnabert.py actually reads.
    trainer_knobs = {"epochs", "batch_size", "encode_batch", "max_length",
                     "learning_rate", "weight_decay"}
    assert set(params) <= trainer_knobs, set(params) - trainer_knobs


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


def test_lone_minority_sample_skips_tuning_instead_of_crashing(features, tmp_path):
    """A single susceptible isolate crashed the tuner instead of skipping.

    With one S among 19 R, both classes exist so the old guard let the search
    run; StratifiedKFold then put the lone S in one validation fold, that
    fold's training half was all-R, and XGBoost raised on single-class y.
    """
    train, test = features
    lone = pd.read_csv(train)
    lone["amikacin"] = "R"
    lone.loc[0, "amikacin"] = "S"
    lopsided = tmp_path / "lopsided.csv"
    lone.to_csv(lopsided, index=False)

    tuned = _tune("xgboost", "amikacin", lopsided, tmp_path / "tune.json")
    assert "skipped" in tuned["tuning"]
    assert tuned["hyperparameters"] == {}

    # Training still fits: one S is enough for a (degenerate) full-data fit.
    metrics = _train("xgboost", "amikacin", lopsided, test, tmp_path / "xgb",
                     params_input=tmp_path / "tune.json")
    assert "skipped" not in metrics


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
