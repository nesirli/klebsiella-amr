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
- `max_samples` default
- reference database locations
- per-tool thread/memory settings

## Makefile targets

| Target | Description |
|--------|-------------|
| `make setup` | Create the two conda environments |
| `make metadata` | Parse metadata and create train/test splits |
| `make qc` | Run fastp + Kraken2 + QUAST |
| `make kraken2` | Run Kraken2 classification |
| `make quast` | Run QUAST assembly statistics |
| `make amr` | Run AMRFinderPlus |
| `make features` | Build gene presence/absence matrix |
| `make sequences` | Build sequence table for DNABERT-2 |
| `make models` | Train XGBoost, LightGBM, and NN |
| `make dnabert` | Train optional DNABERT-2 model |
| `make multiqc` | Build single aggregated MultiQC report |
| `make report` | `multiqc` + aggregated metrics and top features |
| `make all` | `report` |
| `make clean` | Remove generated results, data, and reference databases |

## Disk usage

The pipeline is designed to delete large intermediate files (raw reads, trimmed reads, downsampled reads, assemblies) automatically after downstream steps finish. Small QC and AMR outputs are also removed once the aggregated MultiQC report and feature matrices are built.

- **Reference databases** (~10 GB total) are downloaded once and kept in `reference/`.
- **During a run** you need enough temporary space for reads and assemblies (roughly the size of the input FASTQs plus assembled contigs).
- **After a run** only `results/` is kept: model files, predictions, importance tables, plots, `results/multiqc/multiqc_report.html`, and `results/report/summary.json`.

## Project structure

```text
.
├── Makefile                  # pipeline orchestration
├── config.yaml               # user-editable configuration
├── envs/
│   ├── env-ml.yml            # amr: Python + ML dependencies
│   └── env-bio.yml           # bioinfo: fastp, spades, kraken2, quast, amrfinderplus, multiqc, seqtk
├── metadata.csv              # input metadata (semicolon-delimited)
├── scripts/
│   ├── metadata.py
│   ├── download_reads.py
│   ├── build_features.py
│   ├── build_sequences.py
│   ├── export_config.py
│   ├── summarize.py
│   ├── train_xgboost.py
│   ├── train_lightgbm.py
│   ├── train_nn.py
│   └── train_dnabert.py
├── reference/                # downloaded reference databases
├── data/                     # transient per-sample files (auto-deleted)
└── results/                  # features, models, reports (kept)
```

## Interpretability

Each model emits its own feature-importance CSV and plot:

- **XGBoost / LightGBM**: SHAP summary plots
- **MLP**: permutation importance
- **DNABERT-2**: occlusion importance over gene sequences

The `report` target aggregates per-model metrics and top features across antibiotics.

## Requirements

- [Miniforge](https://github.com/conda-forge/miniforge) (conda + mamba)
- ~10 GB free space for reference databases (Kraken2 standard-8 + AMRFinderPlus)
- Additional temporary space during runs for reads/assemblies

## License

MIT
