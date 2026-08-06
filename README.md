# Klebsiella pneumoniae AMR prediction

A minimal, reproducible pipeline that predicts antimicrobial resistance (AMR) in *Klebsiella pneumoniae* from Illumina reads using tree-based ensemble models and a small neural network.

## Workflow

1. **Metadata** — split NCBI Pathogen Detection metadata by collection year.
2. **Reads** — download paired-end FASTQs from ENA.
3. **QC / trim** — `fastp` trimming and basic reports.
4. **Assembly** — `SPAdes` assembly and custom assembly statistics.
5. **AMR annotation** — `AMRFinderPlus` resistance gene calls.
6. **Features** — gene presence/absence matrix joined to train/test labels.
7. **Models** — XGBoost, LightGBM, and a manual PyTorch MLP, each with per-model interpretability (SHAP or permutation importance).
8. **DNABERT-2** — optional sequence-level model.
9. **Report** — aggregated metrics and top features.

## Quick start

```bash
# 1. Create conda environments
mamba env create -f environment.yml
mamba env create -f environment-bioinfo.yml

# 2. Dev run with 5 samples per split
#    Use 'gmake' on macOS; 'make' on Linux
make metadata MAX_SAMPLES=5
make -j 4 all

# 3. Larger dev run
make metadata MAX_SAMPLES=20
make -j 4 all

# 4. Full run
make metadata MAX_SAMPLES=-1
make -j 8 all
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
| `make reads` | Download reads for all samples |
| `make qc` | Run fastp + MultiQC |
| `make assembly` | Assemble genomes with SPAdes |
| `make amr` | Run AMRFinderPlus |
| `make features` | Build gene presence/absence matrix |
| `make models` | Train XGBoost, LightGBM, and NN |
| `make dnabert` | Train optional DNABERT-2 model |
| `make report` | Aggregate metrics and top features |
| `make all` | `models` + `report` |
| `make clean` | Remove generated results and data |

## Project structure

```text
.
├── Makefile                  # pipeline orchestration
├── config.yaml               # user-editable configuration
├── environment.yml           # amr: Python + ML dependencies
├── environment-bioinfo.yml   # bioinfo: fastp, spades, kraken2, amrfinderplus, multiqc, seqtk
├── metadata.csv              # input metadata (semicolon-delimited)
├── scripts/
│   ├── metadata.py
│   ├── download_reads.py
│   ├── build_features.py
│   ├── build_sequences.py
│   ├── assembly_stats.py
│   ├── export_config.py
│   ├── train_xgboost.py
│   ├── train_lightgbm.py
│   ├── train_nn.py
│   └── train_dnabert.py
├── reference/                # downloaded reference databases
├── data/                     # raw reads, trimmed reads, assemblies
└── results/                  # features, models, reports
```

## Interpretability

Each model emits its own feature-importance CSV and plot:

- **XGBoost / LightGBM**: SHAP summary plots
- **MLP**: permutation importance
- **DNABERT-2**: occlusion importance over gene sequences

The `report` target aggregates per-model metrics and top features across antibiotics.

## Requirements

- [Miniforge](https://github.com/conda-forge/miniforge) (conda + mamba)
- ~8 GB free space for the Kraken2 standard-8 database (downloaded once)
- ~2 GB free space for AMRFinderPlus database (downloaded once)

## License

MIT
