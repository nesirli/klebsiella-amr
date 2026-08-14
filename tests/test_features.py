import json
import subprocess
import sys

import pandas as pd

from test_models import ANTIBIOTICS, _train

GENES_PER_SAMPLE = {
    "SRR00000": ["blaKPC-3", "oqxB"],
    "SRR00001": ["oqxB"],
    "SRR00002": ["blaKPC-3"],
    "SRR00003": ["blaKPC-3", "oqxB"],
    "SRR00004": ["oqxB"],
    "SRR00005": ["blaKPC-3"],
}


def _make_metadata(path, n_train=3, n_test=3):
    rows = []
    for i in range(n_train + n_test):
        year = 2020 if i < n_train else 2023
        rows.append({
            "run": f"SRR{i:05d}",
            "collection_date": f"{year}-01-01",
            "year": year,
            "location": "Test",
            **{abx: "R" if i % 2 == 0 else "S" for abx in ANTIBIOTICS},
        })
    pd.DataFrame(rows).to_csv(path, index=False)


def _make_amr(path, genes):
    pd.DataFrame([{"Element symbol": g} for g in genes]).to_csv(path, sep="\t", index=False)


def _build_features(tmp_path):
    train_meta, test_meta = tmp_path / "train.csv", tmp_path / "test.csv"
    _make_metadata(train_meta)
    _make_metadata(test_meta)

    amr_dir = tmp_path / "amr"
    amr_dir.mkdir()
    for sample, genes in GENES_PER_SAMPLE.items():
        _make_amr(amr_dir / f"{sample}_amr.tsv", genes)

    train_features, test_features = tmp_path / "train_features.csv", tmp_path / "test_features.csv"
    subprocess.run(
        [sys.executable, "scripts/build_features.py",
         "--amr-files", *[str(p) for p in sorted(amr_dir.glob("*.tsv"))],
         "--train-metadata", str(train_meta),
         "--test-metadata", str(test_meta),
         "--train-output", str(train_features),
         "--test-output", str(test_features)],
        check=True,
    )
    return train_features, test_features


def test_build_features_and_xgboost(tmp_path):
    train_features, test_features = _build_features(tmp_path)

    df = pd.read_csv(train_features)
    assert "blaKPC-3" in df.columns
    assert "oqxB" in df.columns

    metrics = _train("xgboost", "amikacin", train_features, test_features, tmp_path / "xgboost")
    assert "f1" in metrics


def test_summarize_collects_every_model_present(tmp_path):
    """dnabert metrics were dropped by a hardcoded model list that omitted it."""
    models_dir = tmp_path / "models"
    for model in ("xgboost", "lightgbm", "nn", "dnabert"):
        model_dir = models_dir / model
        model_dir.mkdir(parents=True)
        (model_dir / "amikacin_metrics.json").write_text(json.dumps({
            "model": model, "antibiotic": "amikacin", "n_train": 10, "n_test": 5,
            "f1": 0.5, "roc_auc": 0.75,
        }))
        pd.DataFrame({"gene": ["blaKPC-3", "oqxB"], "importance": [0.9, 0.1]}).to_csv(
            model_dir / "amikacin_importance.csv", index=False)

    out = tmp_path / "summary.json"
    subprocess.run(
        [sys.executable, "scripts/summarize.py",
         "--models-dir", str(models_dir),
         "--antibiotics", "amikacin",
         "--output", str(out)],
        check=True,
    )

    summary = json.loads(out.read_text())
    assert summary["models"] == ["xgboost", "lightgbm", "nn", "dnabert"]
    assert {row["model"] for row in summary["metrics"]} == {"xgboost", "lightgbm", "nn", "dnabert"}
    assert summary["metrics"][0]["top_features"][0]["gene"] == "blaKPC-3"


def test_summarize_omits_models_that_were_not_run(tmp_path):
    """dnabert is optional; its absence must not fabricate an empty row."""
    models_dir = tmp_path / "models" / "xgboost"
    models_dir.mkdir(parents=True)
    (models_dir / "amikacin_metrics.json").write_text(json.dumps({"f1": 0.5}))

    out = tmp_path / "summary.json"
    subprocess.run(
        [sys.executable, "scripts/summarize.py",
         "--models-dir", str(tmp_path / "models"),
         "--antibiotics", "amikacin",
         "--output", str(out)],
        check=True,
    )

    summary = json.loads(out.read_text())
    assert summary["models"] == ["xgboost"]
    assert len(summary["metrics"]) == 1
