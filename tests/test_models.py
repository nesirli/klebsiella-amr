import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ANTIBIOTICS = ["amikacin", "ciprofloxacin", "ceftazidime", "meropenem"]


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
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)


def _run_model(script, model_output_ext, extra_args, abx, train, test, out_dir):
    out_dir.mkdir(exist_ok=True)
    model_file = out_dir / f"{abx}_model.{model_output_ext}"
    subprocess.run(
        [
            sys.executable,
            script,
            "--train-features",
            str(train),
            "--test-features",
            str(test),
            "--antibiotic",
            abx,
            "--all-antibiotics",
            *ANTIBIOTICS,
            "--model-output",
            str(model_file),
            "--params-output",
            str(out_dir / f"{abx}_params.json"),
            "--metrics-output",
            str(out_dir / f"{abx}_metrics.json"),
            "--predictions-output",
            str(out_dir / f"{abx}_predictions.csv"),
            *extra_args,
        ],
        check=True,
    )
    metrics = json.load(open(out_dir / f"{abx}_metrics.json"))
    for key in ("f1", "precision", "recall", "pr_auc", "confusion_matrix"):
        assert key in metrics
    return metrics


@pytest.mark.parametrize("model", ["xgboost", "lightgbm", "nn"])
def test_model_runs(model):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        train = tmp / "train.csv"
        test = tmp / "test.csv"
        _make_feature_csv(train, n_samples=80)
        _make_feature_csv(test, n_samples=20, seed=1)

        config = {
            "xgboost": ("scripts/train_xgboost.py", "json", [
                "--importance-output", str(tmp / "xgboost" / "{abx}_importance.csv"),
                "--shap-plot-output", str(tmp / "xgboost" / "{abx}_shap.png"),
            ]),
            "lightgbm": ("scripts/train_lightgbm.py", "txt", [
                "--importance-output", str(tmp / "lightgbm" / "{abx}_importance.csv"),
                "--shap-plot-output", str(tmp / "lightgbm" / "{abx}_shap.png"),
            ]),
            "nn": ("scripts/train_nn.py", "pt", [
                "--importance-output", str(tmp / "nn" / "{abx}_importance.csv"),
                "--importance-plot-output", str(tmp / "nn" / "{abx}_importance.png"),
            ]),
        }

        script, ext, extra = config[model]
        out_dir = tmp / model
        for abx in ANTIBIOTICS:
            args = [arg.format(abx=abx) for arg in extra]
            metrics = _run_model(script, ext, args, abx, train, test, out_dir)
            assert metrics["n_train"] == 80
            assert metrics["n_test"] == 20


TUNE_SCRIPTS = {
    "xgboost": ("scripts/tune_xgboost.py", 3),
    "lightgbm": ("scripts/tune_lightgbm.py", 3),
    "nn": ("scripts/tune_nn.py", 2),
}
TRAIN_SCRIPTS = {
    "xgboost": ("scripts/train_xgboost.py", "json", [
        "--importance-output", "{abx}_importance.csv",
        "--shap-plot-output", "{abx}_shap.png",
    ]),
    "lightgbm": ("scripts/train_lightgbm.py", "txt", [
        "--importance-output", "{abx}_importance.csv",
        "--shap-plot-output", "{abx}_shap.png",
    ]),
    "nn": ("scripts/train_nn.py", "pt", [
        "--importance-output", "{abx}_importance.csv",
        "--importance-plot-output", "{abx}_importance.png",
    ]),
}


@pytest.mark.parametrize("model", ["xgboost", "lightgbm", "nn"])
def test_tuning_pipeline(model):
    """End-to-end: tune hyperparameters, then train with tuned params."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        train = tmp / "train.csv"
        test = tmp / "test.csv"
        _make_feature_csv(train, n_samples=80)
        _make_feature_csv(test, n_samples=20, seed=1)

        tune_script, n_trials = TUNE_SCRIPTS[model]
        out_dir = tmp / model
        out_dir.mkdir()

        for abx in ANTIBIOTICS:
            tune_out = out_dir / f"{abx}_tune.json"
            subprocess.run(
                [
                    sys.executable,
                    tune_script,
                    "--train-features",
                    str(train),
                    "--antibiotic",
                    abx,
                    "--all-antibiotics",
                    *ANTIBIOTICS,
                    "--output",
                    str(tune_out),
                    "--n-trials",
                    str(n_trials),
                ],
                check=True,
            )
            tune_result = json.load(open(tune_out))
            assert "best_cv_roc_auc" in tune_result or "skipped" in tune_result

            train_script, ext, extra = TRAIN_SCRIPTS[model]
            model_file = out_dir / f"{abx}_model.{ext}"
            subprocess.run(
                [
                    sys.executable,
                    train_script,
                    "--train-features",
                    str(train),
                    "--test-features",
                    str(test),
                    "--antibiotic",
                    abx,
                    "--all-antibiotics",
                    *ANTIBIOTICS,
                    "--params-input",
                    str(tune_out),
                    "--model-output",
                    str(model_file),
                    "--params-output",
                    str(out_dir / f"{abx}_params.json"),
                    "--metrics-output",
                    str(out_dir / f"{abx}_metrics.json"),
                    "--predictions-output",
                    str(out_dir / f"{abx}_predictions.csv"),
                    *[arg.format(abx=abx) for arg in extra],
                ],
                check=True,
            )
            metrics = json.load(open(out_dir / f"{abx}_metrics.json"))
            assert "f1" in metrics
