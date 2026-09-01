import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, "app")
import predict as P  # noqa: E402

MANIFEST = Path("results/models/manifest.json")
needs_run = pytest.mark.skipif(
    not MANIFEST.exists(), reason="needs artifacts from a completed run")


def test_feature_vector_follows_training_order_not_sample_order():
    """The vector must be indexed by the manifest, never by the upload.

    A model is only meaningful against the column order it was fitted on, and a
    misaligned vector fails silently: it produces a confident number rather
    than an error.
    """
    genes = ["blaKPC-3", "armA", "oqxB"]
    row, present, unknown = P.feature_vector({"oqxB", "blaKPC-3"}, genes)

    assert list(row.columns) == genes
    assert row.iloc[0].tolist() == [1, 0, 1]
    assert present == ["blaKPC-3", "oqxB"]
    assert unknown == []


def test_genes_absent_from_training_are_dropped():
    """A newer AMRFinder database will report genes training never saw; they
    must be reported and discarded, not appended as extra columns."""
    genes = ["armA"]
    row, present, unknown = P.feature_vector({"armA", "brand-new-gene"}, genes)

    assert list(row.columns) == ["armA"]
    assert unknown == ["brand-new-gene"]


def test_rule_fires_on_intrinsic_quinolone_genes():
    """Why the rule is useless for ciprofloxacin (specificity 0.02).

    K. pneumoniae intrinsically carries quinolone efflux genes, so a
    class-presence rule calls nearly every isolate resistant.
    """
    intrinsic = {"NITROFURANTOIN/PHENICOL/QUINOLONE/TIGECYCLINE"}
    assert any(P.RULES["ciprofloxacin"](s) for s in intrinsic)
    # The same isolate implies nothing about carbapenems.
    assert not any(P.RULES["meropenem"](s) for s in intrinsic)


def test_carbapenem_subclass_implies_both_meropenem_and_ceftazidime():
    """Carbapenemases hydrolyse cephalosporins too, so they count for both."""
    assert P.RULES["meropenem"]("CARBAPENEM")
    assert P.RULES["ceftazidime"]("CARBAPENEM")
    assert not P.RULES["amikacin"]("CARBAPENEM")


@needs_run
def test_manifest_columns_match_the_training_matrix():
    manifest = P.load_manifest(MANIFEST)
    features = Path("results/features/train_features.csv")
    if not features.exists():
        pytest.skip("feature matrix removed by 'make report'")
    header = pd.read_csv(features, nrows=0)
    expected = [c for c in header.columns
                if c not in {"run", "collection_date", "year", "location",
                             *manifest["antibiotics"]}]
    assert manifest["genes"] == expected


@needs_run
def test_every_served_model_carries_a_fitted_threshold():
    """Serving at 0.5 would throw away the thresholds; they range 0.31-0.73."""
    manifest = P.load_manifest(MANIFEST)
    assert manifest["models"], "manifest lists no models"
    for per_abx in manifest["models"].values():
        for spec in per_abx.values():
            assert 0.0 < spec["threshold"] < 1.0
            # Paths are relative to the manifest so one file serves both the
            # repo tree and the container image.
            assert not Path(spec["model"]).is_absolute()
            assert (manifest["_base"] / spec["model"]).exists()


@needs_run
def test_serving_reproduces_the_training_prediction():
    """End to end: the app must agree with what training recorded.

    This is the check that catches a feature-alignment regression, which no
    unit test on shapes alone would notice.
    """
    manifest = P.load_manifest(MANIFEST)
    truth = Path("results/models/xgboost/amikacin_predictions.csv")
    if not truth.exists():
        pytest.skip("no trained predictions to compare against")
    rows = pd.read_csv(truth).head(5)

    for _, row in rows.iterrows():
        tsv = Path(f"data/amr/{row['run']}_amr.tsv")
        if not tsv.exists():
            pytest.skip("AMR calls cleaned up")
        symbols, subclasses, _ = P.parse_amr(tsv)
        results, _, _ = P.predict(manifest, symbols, subclasses)
        served = next(e for e in results if e["antibiotic"] == "amikacin")
        assert (served["models"]["xgboost"]["call"] == "R") == bool(row["predicted"])


# Species gate ----------------------------------------------------------------

SKETCH = Path("app/artifacts/kpneumoniae_ref.msh")


def _write_fasta(path, seq, contigs=1):
    with open(path, "w") as fh:
        for i in range(contigs):
            fh.write(f">contig_{i}\n{seq}\n")


def test_reads_and_protein_files_are_rejected_before_anything_expensive(tmp_path):
    """The models only know K. pneumoniae, and AMRFinder costs 30 seconds, so
    an unusable upload has to fail fast and say why."""
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError, match="not a FASTA"):
        P.read_fasta_stats(fastq)

    protein = tmp_path / "prot.faa"
    _write_fasta(protein, "MKVLAWQERTYIPSDFGHKLMKVLAWQERTYIPSDFGHKL")
    with pytest.raises(ValueError, match="does not look like DNA"):
        P.read_fasta_stats(protein)

    empty = tmp_path / "empty.fasta"
    empty.write_text(">contig_only_header\n")
    with pytest.raises(ValueError, match="no sequence"):
        P.read_fasta_stats(empty)


def test_fasta_stats_measure_size_and_gc(tmp_path):
    fasta = tmp_path / "g.fasta"
    _write_fasta(fasta, "GCGC" * 25 + "ATAT" * 25, contigs=4)
    contigs, bases, gc = P.read_fasta_stats(fasta)
    assert contigs == 4
    assert bases == 4 * 200
    assert gc == pytest.approx(0.5)


@pytest.mark.skipif(not SKETCH.exists(), reason="reference sketch not built")
def test_species_gate_rejects_a_non_klebsiella_assembly(tmp_path):
    """Calibration: 60 held-out K. pneumoniae scored at most 0.028 against this
    sketch, E. coli K-12 scored 0.075, and the cutoff sits at 0.05. A random
    sequence must land far outside."""
    import random
    random.seed(0)
    other = tmp_path / "other.fasta"
    _write_fasta(other, "".join(random.choice("ACGT") for _ in range(200_000)))

    accepted, distance, message = P.check_species(other, SKETCH)
    assert not accepted
    assert distance is None or distance > P.SPECIES_MAX_DISTANCE
    assert message


@pytest.mark.skipif(not SKETCH.exists(), reason="reference sketch not built")
def test_species_gate_accepts_a_real_klebsiella_assembly():
    import glob
    assemblies = sorted(glob.glob("data/assembly/*_assembled.fasta"))
    if not assemblies:
        pytest.skip("no assemblies available")

    accepted, distance, _ = P.check_species(assemblies[0], SKETCH)
    assert accepted and distance <= P.SPECIES_MAX_DISTANCE


def test_missing_sketch_fails_closed(tmp_path):
    """A missing sketch must refuse the upload, not silently skip the check --
    otherwise a broken image accepts any organism."""
    fasta = tmp_path / "g.fasta"
    _write_fasta(fasta, "ACGT" * 100)
    accepted, _, message = P.check_species(fasta, tmp_path / "absent.msh")
    assert not accepted and "missing" in message
