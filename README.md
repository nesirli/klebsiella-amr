# Klebsiella pneumoniae AMR prediction

A minimal, reproducible pipeline that predicts antimicrobial resistance (AMR) in *Klebsiella pneumoniae* from Illumina reads using tree-based ensemble models and a small neural network.

## Workflow

1. **Metadata** — split NCBI Pathogen Detection metadata by collection year.
2. **Reads** — download paired-end FASTQs from ENA.
3. **QC / trim** — `fastp` trimming.
4. **Taxonomic QC** — `Kraken2` classification against the standard-8 database.
5. **Assembly QC** — `QUAST` assembly statistics.
6. **Assembly** — `SPAdes` assembly.
7. **AMR annotation** — `AMRFinderPlus` resistance gene calls.
8. **Features** — gene presence/absence matrix joined to train/test labels.
9. **Models** — XGBoost, LightGBM, and a manual PyTorch MLP, each with per-model interpretability (SHAP or permutation importance).
10. **DNABERT-2** — optional sequence-level model.
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

## Makefile targets

| Target | Description |
|--------|-------------|
| `make setup` | Create the two conda environments |
| `make metadata` | Parse metadata and create train/test splits |
| `make tune` | Optuna hyperparameter search (runs automatically before training) |
| `make models` | Train XGBoost, LightGBM, and NN |
| `make dnabert` | Train optional DNABERT-2 model (run before `report` to include it in the summary) |
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

## Requirements

- [Miniforge](https://github.com/conda-forge/miniforge) (conda + mamba)
- ~10 GB free space for reference databases (Kraken2 standard-8 + AMRFinderPlus)
- Additional temporary space during runs for reads/assemblies

## License

MIT
