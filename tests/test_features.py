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
