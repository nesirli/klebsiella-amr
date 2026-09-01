import subprocess
import sys

import pandas as pd
import yaml

ANTIBIOTICS = ["amikacin", "ciprofloxacin"]


def _make_input(path, n_per_year=4):
    """A metadata.csv in NCBI Pathogen Detection shape (semicolon-delimited)."""
    rows = []
    for year in (2019, 2021, 2023, 2024):
        for i in range(n_per_year):
            rows.append({
                "#Run": f"SRR{year}{i:03d}",
                "Collection date": f"{year}-05-01",
                "Location": "Test",
                "AST phenotypes": f"amikacin={'R' if i % 2 else 'S'},"
                                  f"ciprofloxacin={'S' if i % 2 else 'R'}",
            })
    pd.DataFrame(rows).to_csv(path, sep=";", index=False)


def _make_config(path, metadata_path):
    path.write_text(yaml.safe_dump({
        "metadata": str(metadata_path),
        "delimiter": ";",
        "antibiotics": ANTIBIOTICS,
        "splits": {"train_cutoff": 2022, "test_years": [2023, 2024]},
        "max_samples": -1,
        "reference": {"kraken2_db": "reference/kraken2_db",
                      "amrfinder_db": "reference/amrfinderplus"},
        "genome": {"size": 5500000, "target_coverage": 100},
        "resources": {"fastp_threads": 4, "kraken2_threads": 4, "quast_threads": 4,
                      "spades_threads": 8, "spades_memory": 16, "amrfinder_threads": 4},
    }))


def _run_metadata(config, out_dir, max_samples):
    subprocess.run(
        [sys.executable, "scripts/metadata.py",
         "--config", str(config),
         "--train-output", str(out_dir / "train.csv"),
         "--test-output", str(out_dir / "test.csv"),
         "--samples-output", str(out_dir / "samples.txt"),
         "--max-samples", str(max_samples)],
        check=True,
    )
    return (pd.read_csv(out_dir / "train.csv"), pd.read_csv(out_dir / "test.csv"))


def test_temporal_split(tmp_path):
    metadata = tmp_path / "metadata.csv"
    config = tmp_path / "config.yaml"
    _make_input(metadata)
    _make_config(config, metadata)

    train, test = _run_metadata(config, tmp_path, -1)

    assert set(train["year"]) == {2019, 2021}
    assert set(test["year"]) == {2023, 2024}
    assert set(train["run"]).isdisjoint(test["run"])
    assert set(train["amikacin"]) <= {"R", "S"}


def test_max_samples_caps_each_split(tmp_path):
    metadata = tmp_path / "metadata.csv"
    config = tmp_path / "config.yaml"
    _make_input(metadata)
    _make_config(config, metadata)

    train, test = _run_metadata(config, tmp_path, 3)
    assert len(train) == 3
    assert len(test) == 3


def test_max_samples_preserves_class_balance(tmp_path):
    """Capping once drew a near-single-class split from a balanced pool.

    The input is ordered the way NCBI ships it -- all resistant isolates
    first -- so .head(n) returned 20 R and no S. A stratified draw has to
    carry the pool's balance into the capped split.
    """
    metadata = tmp_path / "metadata.csv"
    config = tmp_path / "config.yaml"

    rows = []
    for i in range(80):
        resistant = i < 40  # first half all R: the clustering that broke .head
        rows.append({
            "#Run": f"SRR{i:05d}",
            "Collection date": "2019-05-01",
            "Location": "Test",
            "AST phenotypes": f"amikacin={'R' if resistant else 'S'},"
                              f"ciprofloxacin={'R' if resistant else 'S'}",
        })
    pd.DataFrame(rows).to_csv(metadata, sep=";", index=False)
    _make_config(config, metadata)

    train, _ = _run_metadata(config, tmp_path, 20)

    assert len(train) == 20
    counts = train["amikacin"].value_counts()
    assert counts.get("R", 0) == 10 and counts.get("S", 0) == 10


def test_random_split_puts_both_eras_on_both_sides(tmp_path):
    """The control for the temporal split's distribution shift.

    A temporal split hands the model one era to learn and another to be judged
    on; a stratified random split must mix them, so that a weak score means
    weak features rather than a shifted population.
    """
    metadata = tmp_path / "metadata.csv"
    config = tmp_path / "config.yaml"
    _make_input(metadata, n_per_year=10)
    _make_config(config, metadata)

    subprocess.run(
        [sys.executable, "scripts/metadata.py", "--config", str(config),
         "--train-output", str(tmp_path / "train.csv"),
         "--test-output", str(tmp_path / "test.csv"),
         "--samples-output", str(tmp_path / "samples.txt"),
         "--max-samples", "-1", "--split-mode", "random", "--test-fraction", "0.25"],
        check=True,
    )
    train = pd.read_csv(tmp_path / "train.csv")
    test = pd.read_csv(tmp_path / "test.csv")

    assert set(train["run"]).isdisjoint(test["run"])
    assert len(train) + len(test) == 40
    # Both sides span the full range of years, unlike the temporal split.
    assert set(test["year"]) == {2019, 2021, 2023, 2024}
    assert set(train["year"]) == {2019, 2021, 2023, 2024}
    # And both carry both classes.
    assert set(test["amikacin"]) == {"R", "S"}


def test_unlabelled_samples_are_dropped(tmp_path):
    """No R/S call for any modelled drug means no model can use the isolate."""
    metadata = tmp_path / "metadata.csv"
    config = tmp_path / "config.yaml"

    rows = [{
        "#Run": "SRR00001", "Collection date": "2019-05-01", "Location": "Test",
        "AST phenotypes": "amikacin=R,ciprofloxacin=S",
    }, {
        "#Run": "SRR00002", "Collection date": "2019-05-01", "Location": "Test",
        "AST phenotypes": "meropenem=R",  # not a modelled drug here
    }, {
        "#Run": "SRR00003", "Collection date": "2019-05-01", "Location": "Test",
        "AST phenotypes": "amikacin=ND,ciprofloxacin=ND",
    }]
    pd.DataFrame(rows).to_csv(metadata, sep=";", index=False)
    _make_config(config, metadata)

    train, _ = _run_metadata(config, tmp_path, -1)

    assert list(train["run"]) == ["SRR00001"]
    assert (tmp_path / "samples.txt").read_text().split() == ["SRR00001"]


def test_samples_mk_matches_samples_txt(tmp_path):
    """samples.mk is included by the Makefile; it must list the same runs."""
    metadata = tmp_path / "metadata.csv"
    config = tmp_path / "config.yaml"
    _make_input(metadata)
    _make_config(config, metadata)

    train, test = _run_metadata(config, tmp_path, -1)

    listed = (tmp_path / "samples.txt").read_text().split()
    mk = (tmp_path / "samples.mk").read_text().strip()
    assert mk.startswith("SAMPLES := ")
    assert mk.removeprefix("SAMPLES := ").split() == listed
    assert set(listed) == set(train["run"]) | set(test["run"])


def test_export_config_carries_max_samples(tmp_path):
    """max_samples was documented in config.yaml but never exported to make."""
    metadata = tmp_path / "metadata.csv"
    config = tmp_path / "config.yaml"
    _make_input(metadata)
    _make_config(config, metadata)

    out = tmp_path / "config.mk"
    subprocess.run(
        [sys.executable, "scripts/export_config.py",
         "--config", str(config), "--output", str(out)],
        check=True,
    )

    exported = dict(
        line.split(" := ", 1) for line in out.read_text().splitlines() if " := " in line
    )
    assert exported["MAX_SAMPLES"] == "-1"
    assert exported["KRAKEN_THREADS"] == "4"
    assert exported["ANTIBIOTICS"] == "amikacin ciprofloxacin"
