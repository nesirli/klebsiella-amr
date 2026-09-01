# Klebsiella pneumoniae AMR prediction

A minimal, reproducible pipeline that predicts antimicrobial resistance (AMR) in *Klebsiella pneumoniae* from Illumina reads using tree-based ensemble models and a small neural network.

## Workflow

1. **Metadata** — split NCBI Pathogen Detection metadata by collection year
   (or at random; see *Choosing a split*). Isolates with no R/S call for any
   modelled antibiotic are dropped, since no model can use them.
2. **Reads** — download paired-end FASTQs from ENA.
3. **QC / trim** — `fastp` trimming.
4. **Taxonomic QC** — `Kraken2` classification against the standard-8 database.
5. **Assembly QC** — `QUAST` assembly statistics.
6. **Assembly** — `SPAdes` assembly.
7. **AMR annotation** — `AMRFinderPlus` resistance gene calls.
8. **Features** — gene presence/absence matrix joined to train/test labels.
9. **Models** — XGBoost, LightGBM, and a manual PyTorch MLP, each with per-model interpretability (SHAP or permutation importance).
10. **DNABERT-2** — optional sequence-level model, with its own Optuna tuner.
    Uses the Apple Silicon GPU (MPS) when present, which is ~10x the CPU wheel.
11. **Report** — single aggregated MultiQC report plus model metrics and top features.

## Quick start

```bash
# 1. Create conda environments
mamba env create -f envs/env-ml.yml -n amr
mamba env create -f envs/env-bio.yml -n bioinfo

# 2. Dev run with 5 samples per split (10 samples, 20 FASTQ files)
#    Use 'gmake' on macOS; 'make' on Linux
gmake metadata MAX_SAMPLES=5
gmake -j 4 all

# 3. Larger dev run
gmake metadata MAX_SAMPLES=20
gmake -j 4 all

# 4. Full run
gmake metadata MAX_SAMPLES=-1
gmake -j 8 all
```

## Configuration

Edit `config.yaml` to change:

- input metadata path
- antibiotics and temporal split years
- `max_samples` default (a `MAX_SAMPLES=N` on the command line overrides it)
- reference database locations
- per-tool thread/memory settings

`config.yaml` is exported to make variables by `scripts/export_config.py`, so
changing it re-runs whatever depends on it. Changing `MAX_SAMPLES` re-runs the
split even though it is not a file.

Sizing `resources` to the machine matters more than it looks. The ceiling on
concurrent CPU work is `BATCH_SIZE` × per-tool threads, and `spades_memory` is
a *per-process* cap, so `BATCH_SIZE` × `spades_memory` bounds the assembly
phase. On a 12-core / 18 GB laptop, `BATCH_SIZE=3` with 3 threads and 4 GB per
SPAdes uses about 75% of the cores and stays clear of swap; the original
5 × 8-thread / 16 GB defaults oversubscribed the CPU threefold and exhausted
swap, which was most of the runtime. `model_threads` caps `OMP_NUM_THREADS` for
the model processes — without it XGBoost, LightGBM and torch each grab every
core and fight each other under `make -j`.

## Choosing a split

`splits.mode` decides how train and test are separated, and it changes what a
score means:

- **`temporal`** (default) trains on isolates collected up to `train_cutoff`
  and tests on `test_years`. This is the deployment question — does a model
  fitted on old isolates work on new ones? — and it is the harder, more honest
  number.
- **`random`** holds out a stratified `test_fraction` at random, so both sides
  share a distribution. Use it as a *control*.

Running both is how you tell "the features carry no signal" apart from "the
signal does not survive a shift between eras". On this dataset ceftazidime
scores 0.62–0.68 test ROC AUC under the temporal split and 0.94 under the
random one, from the same features and the same models: its training years are
47% resistant and its test years 85%, because the split is confounded with
which studies deposited isolates when. Amikacin, by contrast, is ~1.0 either
way. A drug that collapses only under the temporal split has a population
problem, not a feature problem.

## Makefile targets

| Target | Description |
|--------|-------------|
| `make setup` | Create the two conda environments |
| `make metadata` | Parse metadata and create train/test splits |
| `make tune` | Optuna hyperparameter search (runs automatically before training) |
| `make models` | Train XGBoost, LightGBM, and NN |
| `make dnabert` | Tune and train the optional DNABERT-2 model (run before `report` to include it in the summary) |
| `make multiqc` | Build single aggregated MultiQC report |
| `make report` | `multiqc` + aggregated metrics; then delete all intermediates |
| `make all` | `report` |
| `make clean` | Remove generated results, data, and reference databases |

## Disk usage

The pipeline deletes intermediate files as soon as they are no longer needed so it can scale to thousands of samples on a laptop.

- **Per-sample cleanup**: as soon as AMRFinderPlus, Kraken2, and QUAST finish for a sample, its raw reads, trimmed reads, and downsampled reads are deleted. Assemblies are kept.
- **Failed samples** are *retired*: every intermediate they own is deleted, including their `fastp` JSON, which drops them out of the active sample list so one bad accession cannot fail the run.
- **`make report` final cleanup**: once the MultiQC report and summary are ready, `make report` deletes everything in `data/` except `data/assembly/`, plus `results/features/` and `results/sequences/`.
- **Reference databases** (~10 GB total) are downloaded once and kept in `reference/`.
- **During a run** you need enough temporary space for the samples currently in flight (roughly the size of their FASTQs plus assembled contigs), multiplied by the `-j` parallelism level.
- **After `make report`** only `results/` and `data/assembly/` are kept: model files, predictions, importance tables, plots, `results/multiqc/multiqc_report.html`, `results/report/summary.json`, and `results/metadata/`.

## Project structure

```text
.
├── Makefile                  # pipeline orchestration
├── config.yaml               # user-editable configuration
├── modules/                  # Makefile includes, one per pipeline stage
│   ├── config.mk             # config.yaml -> make variables
│   ├── metadata.mk           # train/test split
│   ├── pipeline.mk           # per-sample rules: download -> trim -> assemble -> AMR
│   ├── features.mk           # gene matrices
│   ├── models.mk             # tuning + training rules
│   └── batch.mk              # batch driver, MultiQC, report
├── envs/
│   ├── env-ml.yml            # amr: Python + ML dependencies
│   └── env-bio.yml           # bioinfo: fastp, spades, kraken2, quast, amrfinderplus, multiqc, seqtk
├── metadata.csv              # input metadata (semicolon-delimited)
├── scripts/
│   ├── common.py             # shared data loading, metrics, artifacts, Optuna driver
│   ├── mlp.py                # the MLP shared by train_nn.py and tune_nn.py
│   ├── metadata.py
│   ├── download_reads.py
│   ├── build_features.py
│   ├── build_sequences.py
│   ├── export_config.py
│   ├── summarize.py
│   ├── train_{xgboost,lightgbm,nn,dnabert}.py
│   └── tune_{xgboost,lightgbm,nn}.py
├── reference/                # downloaded reference databases
├── data/                     # transient per-sample files (auto-deleted)
└── results/                  # features, models, reports (kept)
```

Every model script shares one contract, implemented in `scripts/common.py`: read a
gene matrix, fit one binary R-vs-S classifier for one antibiotic, and emit the same
six artifacts (model, params, metrics, predictions, importance CSV, importance
plot). Only the estimator and its importance method differ between models. An
antibiotic without enough labelled data to fit still emits all six, marked
`"skipped"`, so a partial dataset cannot fail the run.

Tuners hand hyperparameters to trainers as
`{"hyperparameters": {...}, "tuning": {...}}` — search bookkeeping is kept out of
the block that gets splatted into the estimator.

DNABERT-2 needs its own tuner (`tune_dnabert.py`) rather than the shared
`common.run_tuning`, because it consumes nucleotide strings where the others
consume a float matrix. It reuses the same fold logic and writes the same
handoff file. Its search is deliberately small — 5 trials × 2 folds on a capped
subsample — because every fold refits 117M parameters, so a trial costs minutes
rather than milliseconds.

### Decision thresholds

Each model picks its threshold by maximising balanced accuracy on the
**training** split (`common.choose_threshold`) rather than assuming 0.5, and
records it in the metrics. A fixed 0.5 assumes calibrated probabilities over a
balanced prior, and neither holds here: DNABERT-2 ranked amikacin at 0.980 ROC
AUC yet called all 178 test isolates susceptible, because 74-vs-228 training
labels kept every probability under 0.5.

Caveat worth knowing: the threshold is fitted on predictions the model has
already partly memorised, so it is biased toward the training distribution. For
badly calibrated models this is a large net win; for well-calibrated ones it
can cost a little (XGBoost amikacin balanced accuracy 0.985 → 0.956). Choosing
it on cross-validated out-of-fold predictions instead would remove the bias and
is the obvious next improvement.

## Interpretability

Each model emits a `gene,importance` CSV (most important first) and a plot. The
column is the same across models so they can be compared side by side; what fills
it differs:

- **XGBoost / LightGBM**: mean |SHAP| over the training rows, plotted as a SHAP beeswarm
- **MLP**: permutation importance — drop in test ROC AUC when a gene column is shuffled
- **DNABERT-2**: occlusion importance — drop in test P(resistant) when a gene is removed from the pooled embedding

The `report` target aggregates per-model metrics and the top genes per
model/antibiotic into `results/report/summary.json`. Models that were not run
(DNABERT-2 is optional) are simply absent from it.

## Results from the full run

1184 isolates (970 train / 214 test after the 2020 temporal cutoff), 451 gene
features, four antibiotics, four models. Test-set ROC AUC and balanced
accuracy, with each model's fitted threshold:

| antibiotic | xgboost | lightgbm | nn | dnabert |
|---|---|---|---|---|
| amikacin | 0.997 / 0.96 | 0.997 / 0.96 | **1.000 / 0.99** | 0.947 / 0.82 |
| ciprofloxacin | 0.853 / 0.76 | **0.867 / 0.74** | 0.852 / 0.76 | 0.841 / 0.76 |
| ceftazidime | 0.616 / 0.56 | 0.620 / 0.53 | **0.679 / 0.54** | 0.543 / 0.48 |
| meropenem | 0.866 / 0.82 | **0.882 / 0.86** | 0.867 / 0.82 | 0.781 / 0.71 |

What the numbers say:

- **The biology validates the pipeline.** Amikacin's top genes are
  `aac(6')-Ib`, `armA` and `aph(3')-VIb` — aminoglycoside-modifying enzymes
  and a 16S methyltransferase. Meropenem's are `blaKPC-3`, `blaKPC-2` and
  `ompK36` porin loss. These are the textbook mechanisms, recovered from data.
- **Read `parC_S80I` with suspicion.** It ranks first for all four drugs,
  including ones it has no mechanistic role in. It marks the dominant MDR
  lineage, so it correlates with everything: clonal linkage, not causation.
- **Ceftazidime is a split artifact** (see *Choosing a split*), not a modelling
  failure. Under `mode: random` the same features and models reach 0.94.
- **DNABERT-2 does not pay for itself here.** It matches the tabular models on
  ciprofloxacin and meropenem ranking, trails on amikacin and ceftazidime, and
  costs ~14 hours against seconds for a gradient-boosted tree. Its premise —
  allelic variants a presence flag collapses — is largely already served by
  AMRFinder, which reports point mutations (`gyrA_S83F`, `parC_S80I`) as their
  own columns in the tabular matrix. Its per-run variance is also high: refits
  moved test ROC AUC by up to 0.08.
- **Prefer balanced accuracy to F1 on skewed drugs.** The ceftazidime test
  window is 85% resistant, where calling everything resistant scores F1 ≈ 0.87
  and means nothing.

## Make semantics this pipeline depends on

Non-obvious rules that took a full-scale run to find. Changing them will
silently re-run days of work, or silently delete it.

- **`.SECONDARY:`** (in `Makefile`). GNU make 4.4 deletes files it considers
  intermediate when it exits — here that meant every assembly, AMR table,
  fastp JSON and Kraken report vanishing mid-run, while their grouped-target
  siblings survived. `.SECONDARY` with no prerequisites keeps all of them.
  `.DELETE_ON_ERROR` still removes the target of a recipe that failed.
- **Transient FASTQs are order-only.** Reads, trimmed and downsampled files are
  deleted per sample by design. As *normal* prerequisites, make plans to
  recreate the deleted files, treats that plan as brand-new, and so declares
  every kept artifact stale — sending the analyze phase back to re-download and
  re-assemble every finished sample. After the `|`, they are still built when
  something genuinely needs them. The kept files carry the real timestamps.
- **`download_reads.py` is order-only too.** Editing the downloader must not
  invalidate reads already on disk. To genuinely refetch a sample, delete its
  reads.
- **MultiQC depends on the `cleaned` markers, not on `qc/.done`.** Reaching for
  a fastp JSON drags in the grouped rule that produced it, and that group also
  lists the trimmed FASTQs which cleanup has deleted — an incomplete group, so
  make re-runs fastp for every sample. A `cleaned` marker proves QC, Kraken2
  and QUAST all finished and belongs to no group.
- **A capped run must sample, not truncate.** `MAX_SAMPLES` uses a stratified
  draw; NCBI metadata arrives grouped by submitting study, so `.head(n)` once
  produced a split that was 18 R and 1 S for a drug whose pool is 51% R.

## Failure handling

Transient problems must not be mistaken for bad data, because the pipeline's
response to bad data is permanent: a failed sample is *retired* and drops out
of the run.

- `download_reads.py` separates "cannot reach ENA" (DNS failures, resets,
  timeouts, HTTP 403/429/5xx) from "this accession is bad" (404, malformed
  filereport). Network failures do not consume the retry budget; they back off
  and wait up to 30 minutes. A brief outage once retired 57 healthy samples.
- A truncated-but-advancing transfer counts as progress, not failure. ENA
  closes long transfers early, and one 186 MB file needed several resumed
  passes; charging each pass against the retry budget retired a sample that was
  downloading correctly.
- Downloads are verified by decompressing them, not by size. One FASTQ arrived
  at exactly its advertised `Content-Length` and still ended mid-stream — only
  fastp noticed, three steps later, and the sample was retired for it.
- Retired samples are listed in `results/metadata/retired.txt`. **Re-running
  the same command recovers them**, because finished samples are skipped via
  their `data/cleaned/*_cleaned` markers. Check that file when a run ends.

## Known issues

- `build_sequences.py` keys genes on the third whitespace field of the
  AMRFinder FASTA header, which its docstring calls "the SAME key
  build_features.py uses". It is not: 378 of 686 keys are coordinate-stamped
  frameshift entries such as `cirA_@del_240_241_24418_24421_0_STOP`, unique per
  isolate, where the tabular matrix has 451 clean `Element symbol` values. So
  DNABERT-2 sees hundreds of near-singleton pseudo-genes. Mapping the header
  back to the TSV's `Element symbol` would fix it, and would require re-tuning.
- Kraken2 runs for every sample and lands in the MultiQC report, but nothing
  consumes it. Contaminated or misidentified isolates are not excluded.
- `make report` deletes `data/amr/`, so adding DNABERT-2 afterwards means
  re-running AMRFinder over the kept assemblies (a couple of hours, no
  re-assembly). Run `make dnabert` before `make report` to avoid that.

## Serving

`app/` is a self-contained Streamlit service: upload an assembled genome, get a
resistance profile. `make app-artifacts` then `make app`, or build the
Dockerfile in that folder. `app/README.md` covers deployment; the notes below
are about why it is built the way it is.

**Why it shows a rule and a model side by side.** Neither wins everywhere.
Measured on the same test set, balanced accuracy:

| antibiotic | rule | xgboost | lightgbm |
|---|---|---|---|
| amikacin | 0.703 | 0.956 | 0.959 |
| ciprofloxacin | 0.508 | 0.761 | 0.741 |
| ceftazidime | 0.509 | 0.561 | 0.531 |
| meropenem | **0.928** | 0.821 | 0.856 |

The rule is "resistant if AMRFinder reports any gene of that drug class". It
wins for meropenem, because carbapenemase carriage is a clean determinant. It
is useless for ciprofloxacin and ceftazidime — specificity 0.02 and 0.04 —
because *K. pneumoniae* intrinsically carries quinolone efflux genes
(`oqxA`/`oqxB`) and a chromosomal beta-lactamase, so "has a gene of this class"
is true of nearly every isolate. The models learn to ignore those and weigh the
QRDR mutations that matter. AMRFinder runs either way: it is the models' input,
not an alternative to them.

**The manifest is the correctness boundary.** `results/models/manifest.json`
freezes the gene columns *in training order*, plus each model's fitted
threshold. Serving must index features by the manifest, never by the uploaded
sample — a misaligned vector produces a confident number rather than an error.
Genes the sample has that training lacked are dropped and reported; a large
count there means the AMRFinder database has moved since training.
`test_predict.py` asserts a served prediction reproduces the training
prediction, which is what catches an alignment regression.

Model paths inside the manifest are stored **relative to the manifest file**,
so the same document serves the repo tree (`results/models/`) and the image
(`/app/artifacts/`) with no rewriting; `predict.py` resolves them against the
manifest's own directory.

Build the manifest **before** `make report`: it needs the feature matrices,
which that target deletes.

The image cannot use the `python:slim` + `uv` pattern of the other portfolio
services, because AMRFinderPlus is a bioconda package that shells out to BLAST+
and HMMER. micromamba installs both stacks from one spec. `predict.py` calls
`amrfinder` directly when it is on PATH (the container) and falls back to
`conda run -n bioinfo` otherwise (a local checkout).

**Species gate.** The models were fitted only on *K. pneumoniae*, so the app
refuses anything else rather than scoring it. Three checks run before
AMRFinder, cheapest first:

1. Is it nucleotide FASTA at all? Catches uploaded FASTQ reads and protein
   FASTA, which would otherwise fail obscurely several tools later.
2. Is the assembly 3-9 Mb? Catches single contigs, read files and eukaryotic
   sequence without running anything expensive.
3. Mash distance to `artifacts/kpneumoniae_ref.msh`, a 1000-hash sketch of 20
   assemblies from this dataset (164 KB).

The 0.05 cutoff was calibrated, not guessed: 60 held-out *K. pneumoniae*
assemblies scored at most **0.028** (median 0.009), while *E. coli* K-12 --
the hardest realistic confusion, same family -- scored **0.075**. The cutoff
sits in that gap and is close to the conventional 95% ANI species boundary.
Members of the *K. pneumoniae* complex (*K. variicola*, *K. quasipneumoniae*)
fall under it and are accepted, which is intended: they share the AMR genetics.

It **fails closed**. A missing sketch or a missing `mash` binary rejects the
upload rather than skipping the check, so a broken image cannot silently start
accepting any organism. `docker-entrypoint.sh` also asserts the sketch exists
at start-up.

This is the gate the pipeline itself still lacks -- Kraken2 runs over every
training sample but nothing consumes its verdict, so a contaminated isolate can
still reach the training set.

**Not served:** the MLP (adds little over the trees) and DNABERT-2 (450 MB per
antibiotic, slow, and worse on every drug).

**The password is a gate, not authentication** — one shared secret, no
accounts, no lockout, no audit trail. It is compared with `hmac.compare_digest`
and read from `.streamlit/secrets.toml` (gitignored) or `$AMR_APP_PASSWORD`.
Serve over HTTPS or it crosses the wire in clear text, and put a real proxy in
front before widening access. Uploads go to a temporary directory that is
deleted when the analysis finishes, and are never logged.

## Requirements

- [Miniforge](https://github.com/conda-forge/miniforge) (conda + mamba)
- ~10 GB free space for reference databases (Kraken2 standard-8 + AMRFinderPlus)
- Additional temporary space during runs for reads/assemblies

## License

MIT
