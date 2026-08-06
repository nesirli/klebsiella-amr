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
    assert "f1" in metrics
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
