import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from test_models import ANTIBIOTICS


def _make_metadata(path, n_train=3, n_test=3):
    rows = []
    for i in range(n_train):
        rows.append({
            "run": f"SRR{i:05d}",
            "collection_date": "2020-01-01",
            "year": 2020,
            "location": "Test",
            **{abx: "R" if i % 2 == 0 else "S" for abx in ANTIBIOTICS},
        })
    for i in range(n_test):
        rows.append({
            "run": f"SRR{i + n_train:05d}",
            "collection_date": "2023-01-01",
            "year": 2023,
            "location": "Test",
            **{abx: "R" if i % 2 == 0 else "S" for abx in ANTIBIOTICS},
        })
    pd.DataFrame(rows).to_csv(path, index=False)


def _make_amr(path, genes):
    rows = [{"Element symbol": g} for g in genes]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_build_features_and_xgboost():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        train_meta = tmp / "train.csv"
        test_meta = tmp / "test.csv"
        _make_metadata(train_meta)
        _make_metadata(test_meta)

        amr_dir = tmp / "amr"
        amr_dir.mkdir()
        genes_per_sample = {
            "SRR00000": ["blaKPC-3", "oqxB"],
            "SRR00001": ["oqxB"],
            "SRR00002": ["blaKPC-3"],
            "SRR00003": ["blaKPC-3", "oqxB"],
            "SRR00004": ["oqxB"],
            "SRR00005": ["blaKPC-3"],
        }
        for sample, genes in genes_per_sample.items():
            _make_amr(amr_dir / f"{sample}_amr.tsv", genes)

        train_features = tmp / "train_features.csv"
        test_features = tmp / "test_features.csv"
        subprocess.run(
            [
                sys.executable,
                "scripts/build_features.py",
                "--amr-files",
                *[str(p) for p in amr_dir.glob("*.tsv")],
                "--train-metadata",
                str(train_meta),
                "--test-metadata",
                str(test_meta),
                "--train-output",
                str(train_features),
                "--test-output",
                str(test_features),
            ],
            check=True,
        )

        df = pd.read_csv(train_features)
        assert "blaKPC-3" in df.columns
        assert "oqxB" in df.columns

        out_dir = tmp / "xgboost"
        out_dir.mkdir()
        subprocess.run(
            [
                sys.executable,
                "scripts/train_xgboost.py",
                "--train-features",
                str(train_features),
                "--test-features",
                str(test_features),
                "--antibiotic",
                "amikacin",
                "--all-antibiotics",
                *ANTIBIOTICS,
                "--model-output",
                str(out_dir / "model.json"),
                "--params-output",
                str(out_dir / "params.json"),
                "--metrics-output",
                str(out_dir / "metrics.json"),
                "--predictions-output",
                str(out_dir / "predictions.csv"),
                "--importance-output",
                str(out_dir / "importance.csv"),
                "--shap-plot-output",
                str(out_dir / "shap.png"),
            ],
            check=True,
        )

        metrics = json.load(open(out_dir / "metrics.json"))
        assert "f1" in metrics
