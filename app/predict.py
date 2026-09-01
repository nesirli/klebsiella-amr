#!/usr/bin/env python3
"""Predict a resistance profile for one assembled genome.

Two independent answers per antibiotic, deliberately shown side by side:

- **rule**: does AMRFinder report any gene of that drug class? Transparent, and
  on this dataset it beats the models for meropenem (balanced accuracy 0.93 vs
  0.86), because carbapenemase carriage is a clean determinant.
- **model**: XGBoost / LightGBM over the gene presence matrix. Much better
  where the rule collapses -- K. pneumoniae intrinsically carries quinolone
  efflux genes (oqxA/oqxB) and a chromosomal beta-lactamase, so "has a
  quinolone gene" is true of nearly every isolate and the rule calls everything
  resistant (specificity 0.02 for ciprofloxacin).

Where the two disagree is where a human should look, which is the reason for
showing both rather than picking one.

Feature alignment comes from the manifest, never from the uploaded sample: the
vector must have exactly the training columns, in training order.
"""
import argparse
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

# Species gate. The models were fitted only on K. pneumoniae, so anything else
# would get a confident answer drawn from the wrong biology.
#
# Threshold calibrated against this dataset: 60 held-out K. pneumoniae
# assemblies scored a maximum Mash distance of 0.028 to the reference sketch
# (median 0.009), while E. coli K-12 -- the hardest realistic confusion, same
# family -- scored 0.075. 0.05 sits in the gap with room on both sides, and
# corresponds to roughly the conventional 95% ANI species boundary.
SPECIES_SKETCH = "kpneumoniae_ref.msh"
SPECIES_MAX_DISTANCE = 0.05

# A Klebsiella assembly is ~5.0-5.9 Mb. This wider band only has to reject
# things that are not a bacterial genome at all -- read files, a single contig,
# a eukaryotic chromosome -- before spending 30 seconds in AMRFinder.
MIN_GENOME_BP = 3_000_000
MAX_GENOME_BP = 9_000_000

# A hit implicates the drug when its AMRFinder Subclass matches. Substring
# matching on purpose: subclasses are slash-joined lists like
# "AMIKACIN/KANAMYCIN/TOBRAMYCIN".
RULES = {
    "amikacin": lambda s: "AMIKACIN" in s,
    "ciprofloxacin": lambda s: "QUINOLONE" in s,
    "ceftazidime": lambda s: "CEPHALOSPORIN" in s or "CARBAPENEM" in s,
    "meropenem": lambda s: "CARBAPENEM" in s,
}

# Drugs whose temporal-split score says the model does not transfer to new
# isolates, even though it scores ~0.94 within the training distribution.
LOW_CONFIDENCE = {"ceftazidime"}


def load_manifest(path):
    """Model paths inside the manifest are relative to the manifest itself, so
    the same file serves the repo tree and the image without edits."""
    manifest = json.loads(Path(path).read_text())
    manifest["_base"] = Path(path).parent
    return manifest


def tool_command(tool):
    """How to invoke a bioinformatics tool here.

    In the container these are installed into the base environment and are on
    PATH. In a local checkout they live in the `bioinfo` conda environment
    while the app runs from `amr`, so they need `conda run`.
    """
    if shutil.which(tool):
        return [tool]
    env = os.environ.get("BIOINFO_ENV", "bioinfo")
    return ["conda", "run", "--no-capture-output", "-n", env, tool]


def amrfinder_command():
    return tool_command("amrfinder")


def read_fasta_stats(path):
    """-> (contigs, total_bp, gc_fraction). Raises ValueError if not nucleotide
    FASTA, which is the common case of someone uploading reads or a protein
    file and getting an inscrutable failure three tools later."""
    contigs, bases, gc, acgtn = 0, 0, 0, 0
    with open(path, errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                contigs += 1
                continue
            if contigs == 0:
                raise ValueError("not a FASTA file (no '>' header before sequence)")
            upper = line.upper()
            bases += len(upper)
            gc += upper.count("G") + upper.count("C")
            acgtn += sum(upper.count(c) for c in "ACGTN")
    if contigs == 0:
        raise ValueError("not a FASTA file (no '>' header found)")
    if bases == 0:
        raise ValueError("FASTA contains headers but no sequence")
    if acgtn / bases < 0.9:
        raise ValueError("sequence does not look like DNA (is this a protein FASTA?)")
    return contigs, bases, gc / bases


def check_species(fasta, sketch, max_distance=SPECIES_MAX_DISTANCE):
    """Mash distance from the assembly to a K. pneumoniae reference sketch.

    -> (accepted, distance, explanation). A missing sketch or a missing mash
    binary is reported rather than silently skipped: quietly disabling the gate
    would let the wrong organism through with a confident answer.
    """
    if not Path(sketch).exists():
        return False, None, f"species reference sketch missing at {sketch}"
    try:
        proc = subprocess.run(
            [*tool_command("mash"), "dist", str(sketch), str(fasta)],
            check=True, capture_output=True, text=True,
        )
    except FileNotFoundError:
        return False, None, "mash is not installed, so species cannot be verified"
    except subprocess.CalledProcessError as exc:
        return False, None, f"mash failed: {exc.stderr.strip()[:200]}"

    distances = []
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            try:
                distances.append(float(parts[2]))
            except ValueError:
                continue
    if not distances:
        return False, None, "mash returned no comparisons"

    best = min(distances)
    if best <= max_distance:
        return True, best, f"Mash distance {best:.3f} to K. pneumoniae reference"
    return False, best, (
        f"Mash distance {best:.3f} exceeds the {max_distance} species cutoff "
        f"(K. pneumoniae assemblies score below 0.03; E. coli scores 0.075)")


def run_amrfinder(fasta, workdir, database, threads=4):
    """Run AMRFinderPlus on an assembly; returns the TSV path."""
    out = Path(workdir) / "amr.tsv"
    subprocess.run(
        [*amrfinder_command(),
         "--nucleotide", str(fasta),
         "--organism", "Klebsiella_pneumoniae",
         "--database", str(Path(database) / "latest"),
         "--output", str(out),
         "--threads", str(threads)],
        check=True, capture_output=True, text=True,
    )
    return out


def parse_amr(tsv):
    """-> (gene symbols, subclass strings, table for display)."""
    df = pd.read_csv(tsv, sep="\t")
    symbols = set(df["Element symbol"].dropna().astype(str))
    subclasses = set(df["Subclass"].dropna().astype(str))
    columns = [c for c in ("Element symbol", "Element name", "Class", "Subclass",
                           "% Identity to reference") if c in df.columns]
    return symbols, subclasses, df[columns]


def feature_vector(symbols, genes):
    """One row, exactly the manifest's columns in the manifest's order.

    Genes the sample has that training never saw are dropped; genes training
    had that the sample lacks are zero. Both are expected, and the count of
    dropped genes is worth surfacing: a large number means the AMRFinder
    database has moved since training.
    """
    present = [g for g in genes if g in symbols]
    row = pd.DataFrame([[1 if g in symbols else 0 for g in genes]], columns=genes)
    return row, present, sorted(symbols - set(genes))


def _load_booster(name, path):
    if name == "xgboost":
        import xgboost
        booster = xgboost.Booster()
        booster.load_model(path)
        return lambda x: booster.predict(xgboost.DMatrix(x))
    if name == "lightgbm":
        import lightgbm
        booster = lightgbm.Booster(model_file=path)
        return lambda x: booster.predict(x)
    raise ValueError(f"unsupported model: {name}")


def predict(manifest, symbols, subclasses):
    """Rule and model calls per antibiotic."""
    row, present, unknown = feature_vector(symbols, manifest["genes"])
    results = []
    for abx in manifest["antibiotics"]:
        rule = RULES[abx]
        rule_hits = sorted(s for s in subclasses if rule(s))
        entry = {
            "antibiotic": abx,
            "rule_call": "R" if rule_hits else "S",
            "rule_evidence": rule_hits,
            "models": {},
            "low_confidence": abx in LOW_CONFIDENCE,
        }
        for name, per_abx in manifest["models"].items():
            spec = per_abx.get(abx)
            if not spec:
                continue
            model_path = Path(spec["model"])
            if not model_path.is_absolute():
                model_path = Path(manifest.get("_base", ".")) / model_path
            proba = float(np.asarray(
                _load_booster(name, str(model_path))(row)).ravel()[0])
            entry["models"][name] = {
                "probability": proba,
                "call": "R" if proba >= spec["threshold"] else "S",
                "threshold": spec["threshold"],
            }
        results.append(entry)
    return results, present, unknown


def parse_args():
    p = argparse.ArgumentParser(description="Predict resistance for one genome")
    p.add_argument("--fasta", help="assembled genome (skipped if --amr-tsv given)")
    p.add_argument("--amr-tsv", help="precomputed AMRFinder output")
    p.add_argument("--manifest", required=True)
    p.add_argument("--amrfinder-db", default="reference/amrfinderplus")
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args()


def main():
    args = parse_args()
    manifest = load_manifest(args.manifest)

    if args.amr_tsv:
        symbols, subclasses, _ = parse_amr(args.amr_tsv)
    elif args.fasta:
        with tempfile.TemporaryDirectory() as tmp:
            tsv = run_amrfinder(args.fasta, tmp, args.amrfinder_db, args.threads)
            symbols, subclasses, _ = parse_amr(tsv)
    else:
        raise SystemExit("need --fasta or --amr-tsv")

    results, present, unknown = predict(manifest, symbols, subclasses)
    print(json.dumps({"genes_matched": len(present),
                      "genes_not_in_training": len(unknown),
                      "predictions": results}, indent=2))


if __name__ == "__main__":
    main()
