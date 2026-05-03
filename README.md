# Interpretable Deep-Learning and Ensemble Models for Predicting Multidrug Resistance in *Klebsiella pneumoniae*

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.32-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Python](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/)

A comprehensive, reproducible Snakemake workflow for genomic prediction of antimicrobial resistance (AMR) in *Klebsiella pneumoniae* using tree-based ensemble methods and deep learning architectures with temporal validation and interpretability analysis.

## 📋 Overview

This pipeline implements a 20-stage workflow comparing four machine learning architectures (XGBoost, LightGBM, 1D CNN, DNABERT-2) for predicting resistance to four critical antibiotic classes:
- **Carbapenems** (meropenem)
- **Cephalosporins** (ceftazidime)
- **Fluoroquinolones** (ciprofloxacin)
- **Aminoglycosides** (amikacin)

**Key Features:**
- ✅ Rigorous temporal validation (pre-2023 training → 2023-2024 testing)
- ✅ SHAP-based interpretability analysis
- ✅ Comprehensive quality control pipeline
- ✅ Automated feature engineering and selection
- ✅ Full reproducibility with conda environments
- ✅ Optimized for high-performance computing


## 🚀 Quick Start

### Prerequisites

**Hardware Requirements:**
- **Recommended:** 32 vCPUs, 128GB RAM, 1TB SSD (Hetzner server or equivalent)
- **Minimum (Mac/laptop):** 8 cores, 16GB RAM (reduced parallelism, longer runtime)
- **Runtime:** 2-3 days for complete pipeline on recommended hardware
- **Storage:** ~700GB for metadata and results

**Software (choose one):**
- **Option A — Docker** (recommended, works on any OS including Mac ARM64):
  - [Docker Desktop](https://www.docker.com/products/docker-desktop/) ≥24.0
  - [Docker Compose](https://docs.docker.com/compose/) V2
- **Option B — Native Conda:**
  - [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/) (recommended)
  - [Snakemake](https://snakemake.readthedocs.io/) ≥7.32
  - Git

### Option A: Docker (Recommended)

Docker eliminates all dependency/platform issues — bioinformatics tools like SPAdes, freebayes, and kraken2 run in a Linux x86_64 container regardless of host OS.

```bash
# Clone repository
git clone https://github.com/nesirli/msc-project.git
cd msc-project

# Build the image (first time takes ~30-60 min to create all conda envs)
docker compose build

# Run the full pipeline
docker compose up pipeline

# Or run individual stages
docker compose up preprocess     # Stages 1-5: QC + assembly
docker compose up train          # Stages 14-18: all models

# Interactive shell inside the container
docker compose run --rm dev

# Run tests
docker compose run --rm test
```

**Adjusting resources (edit `.env` or pass as environment):**
```bash
# Use 4 threads and limit to 32 GB RAM
THREADS=4 DOCKER_MEMORY=32G docker compose up pipeline

# On Mac with 8 cores:
THREADS=8 DOCKER_CPUS=8 DOCKER_MEMORY=16G docker compose up pipeline
```

**GPU training (NVIDIA only):**
```bash
# Requires nvidia-container-toolkit
docker compose --profile gpu up train-gpu
```

### Option B: Native Installation (Conda)

```bash
# Clone repository
git clone https://github.com/nesirli/msc-project.git
cd msc-project

# Install Snakemake (if not already installed)
conda create -n snakemake -c conda-forge -c bioconda snakemake=7.32
conda activate snakemake
```

### Running the Pipeline

**Full Pipeline (20 stages):**
```bash
# Maximum parallelization strategy (recommended for 32 vCPU)
bash run_max_parallel.sh

# OR standard Snakemake execution
snakemake --use-conda --cores 32 --jobs 16
```

**Individual Stages:**
```bash
# Run specific stage independently
snakemake --use-conda --cores 32 -s rules/01_metadata.smk metadata_all
snakemake --use-conda --cores 32 -s rules/06_amr_analysis.smk amr_analysis_all
```

**Partial Workflows:**
```bash
# Preprocessing only (stages 1-5)
snakemake --use-conda --cores 32 preprocess

# Feature extraction (stages 6-10)
snakemake --use-conda --cores 32 feature_extraction

# Tree models only (stages 14-15)
snakemake --use-conda --cores 32 tree_models

# Deep learning models (stages 16-18)
snakemake --use-conda --cores 32 dl_models

# Interpretability analysis (stage 19)
snakemake --use-conda --cores 32 interpretability
```

### Running Tests

```bash
# Install test dependencies
python -m venv venv
source venv/bin/activate
pip install -r requirements-dev.txt

# Run all tests
python -m pytest tests/ -v

# Run with coverage report
python -m pytest tests/ --cov=utils --cov-report=term-missing
```


## 📊 Pipeline Stages

The workflow consists of 20 interconnected stages, organized into 5 functional groups:

### Stage 1-5: Data Acquisition & Quality Control

**Stage 1: Metadata Processing** (`rules/01_metadata.smk`)
- Parse NCBI pathogen metadata CSV
- Validate phenotype data and filter out incomplete records
- Generate train/test splits using temporal stratification (pre-2023 vs 2023-2024)
- Output: `metadata_train_processed.csv`, `metadata_test_processed.csv`

**Stage 2: Download** (`rules/02_download.smk`)
- Retrieve paired-end Illumina reads from SRA using fastq-dump
- Verify integrity with md5 checksums
- Parallel download with configurable concurrency
- Output: `data/reads/{run_accession}_{1,2}.fastq.gz`

**Stage 3: Pre-assembly QC** (`rules/03_preassembly_qc.smk`)
- FastQC quality analysis
- fastp trimming (adapters, low-quality bases, short reads <50bp)
- Contamination screening with Kraken2
- Output: `results/qc/preassembly_multiqc.html`

**Stage 4: Assembly** (`rules/04_assembly.smk`)
- SPAdes de novo genome assembly (k-mer sizes: auto)
- Configurable kmer selection and memory usage
- Output: `assemblies/{run_accession}/contigs.fasta`

**Stage 5: Post-assembly QC** (`rules/05_postassembly_qc.smk`)
- QUAST assembly metrics (N50, GC%, contig count)
- Kraken2 species/contamination confirmation
- MultiQC HTML report aggregation
- Output: `results/qc/postassembly_multiqc.html`

### Stage 6-10: Feature Engineering

**Stage 6: AMR Analysis** (`rules/06_amr_analysis.smk`)
- AMRFinderPlus resistance gene annotation
- Extract presence/absence of 12,000+ known resistance genes
- Database-aware calling with curated gene families
- Output: `results/amr/combined_amrfinder.csv`

**Stage 7: SNP Analysis** (`rules/07_snp_analysis.smk`)
- Map reads to reference genome (BWA)
- Call SNPs with Snippy (min quality 30, min depth 10)
- Extract core-genome variants ~1.2M SNPs
- Output: `results/snp/core_genome_variants.vcf`

**Stage 8: Feature Matrix** (`rules/08_feature_matrix.smk`)
- Combine AMR genes + SNPs into unified matrix
- Encode as binary (0/1) features
- 1.2M total features × 1,372 samples
- Output: `results/features/feature_matrix_raw.csv`

**Stage 9: Feature Selection** (`rules/09_feature_selection.smk`)
- Chi-square test for association with phenotype
- Mutual Information scoring
- Select top 325 features (99.97% dimensionality reduction)
- Output: `results/features/{antibiotic}_feature_importance.csv`

**Stage 10: Batch Correction** (`rules/10_batch_correction.smk`)
- Detect geographic/temporal batch effects
- ComBat-Seq empirical Bayes correction
- Preserve biological signal while removing batch
- Output: `results/batch_correction/features_batch_corrected.csv`

### Stage 11-13: Dataset Preparation for ML

**Stage 11: Training Data Preparation** (`rules/11_prepare_training_data.smk`)
- Split into train (pre-2023) and test (2023-2024) sets
- Stratify by geographic origin and phenotype
- SMOTE-ENN resampling on training set only
- Output: `results/features/tree_models/{antibiotic}_{train,test}_final.csv`

**Stage 12: K-mer Datasets** (`rules/12_create_kmer_datasets.smk`)
- Tokenize DNA sequences into k-mers (k=11)
- Create one-hot encoded input for 1D CNN
- Store as binary NPZ format (space-efficient)
- Output: `results/features/deep_models/{antibiotic}_kmer_{train,test}_final.npz`

**Stage 13: DNABERT Datasets** (`rules/13_create_dnabert_datasets.smk`)
- Tokenize DNA sequences for DNABERT-2 tokenizer
- Create input_ids and attention_mask tensors
- Prepare for Hugging Face transformers pipeline
- Output: `results/features/deep_models/{antibiotic}_dnabert_{train,test}_final.npz`

### Stage 14-18: Model Training & Evaluation

**Stage 14: XGBoost** (`rules/14_train_xgboost.smk`)
- Gradient boosting tree model on selected features
- GridSearchCV hyperparameter optimization (5-fold CV)
- SHAP feature importance computation
- Output: `results/models/xgboost/{antibiotic}_results.json`

**Stage 15: LightGBM** (`rules/15_train_lightgbm.smk`)
- Light gradient boosting alternative to XGBoost
- Faster training, often similar performance
- Nested CV with hyperparameter tuning
- Output: `results/models/lightgbm/{antibiotic}_results.json`

**Stage 16: 1D CNN** (`rules/16_train_1dcnn.smk`)
- Convolutional neural network on k-mer spectra
- Architecture: Conv1D → MaxPool → Dense → Dropout
- Adam optimizer, binary crossentropy loss
- Output: `results/models/cnn/{antibiotic}_results.json`

**Stage 17: Sequence CNN** (`rules/17_train_sequence_cnn.smk`)
- CNN on raw DNA sequences (ACGT encoding)
- Larger receptive field, fewer parameters than 1D CNN
- Early stopping on validation set
- Output: `results/models/sequence_cnn/{antibiotic}_results.json`

**Stage 18: DNABERT** (`rules/18_train_dnabert.smk`)
- Fine-tune pre-trained DNA transformer (DNABERT-2)
- Sequence length 512 bp with sliding window approach
- Hugging Face trainer with linear warmup schedule
- Output: `results/models/dnabert/{antibiotic}_results.json`

### Stage 19-20: Analysis & Interpretation

**Stage 19: Interpretability Analysis** (`rules/19_interpretability_analysis.smk`)
- SHAP TreeExplainer for tree models
- SHAP DeepExplainer for neural networks
- Generate feature importance plots (bar, dependence, beeswarm)
- Biological validation using motif databases
- Output: `results/interpretability/{antibiotic}_*.png` and `*.json`

**Stage 20: Ensemble Analysis** (`rules/20_ensemble_analysis.smk`)
- Load predictions from all 5 models
- Evaluate ensemble strategies (average, weighted, voting)
- Compare ensemble vs best single model
- Generate final performance comparison
- Output: `results/ensemble/{antibiotic}_ensemble_analysis.json`

## 📁 Project Structure

```
msc-project/
├── Snakefile                    # Master workflow orchestrator (20 stages)
├── config/
│   └── config.yaml              # Pipeline configuration
├── rules/                       # Individual Snakemake rule modules
│   ├── 01_metadata.smk          # Metadata curation
│   ├── 02_download.smk          # SRA data retrieval
│   ├── 06_amr_analysis.smk      # Resistance gene annotation
│   ├── 14_train_xgboost.smk     # XGBoost model training
│   └── ...20 total stages
├── scripts/                     # Python implementation scripts
│   ├── 01_metadata.py
│   ├── 14_train_xgboost.py
│   ├── 20_ensemble_analysis.py  # Ensemble evaluation
│   └── 21_ensemble_summary.py   # Summary statistics
├── envs/                        # Conda environment specifications (16 environments)
│   ├── assembly.yaml            # SPAdes, QUAST, Kraken2
│   ├── xgboost.yaml             # XGBoost + SHAP
│   ├── dnabert.yaml             # Hugging Face transformers
│   └── ...
├── utils/                       # Shared utility modules
│   ├── cross_validation.py      # Temporal and geographic CV strategies
│   ├── evaluation.py            # Comprehensive metrics computation
│   ├── ensemble_methods.py      # Model combination strategies
│   ├── motif_analysis.py        # Biological feature validation
│   └── output_validation.py     # Result standardization
├── tests/                       # Unit and integration tests
│   ├── test_dl_training.py
│   ├── test_ensemble_methods.py
│   └── ...
├── data/
│   ├── metadata.csv             # User-provided NCBI pathogen metadata
│   └── reference/               # Reference genome and databases
│       ├── reference_genome     # K. pneumoniae HS11286
│       ├── genome_size.txt
│       └── kraken2_db/          # Contamination screening database
├── results/                     # Auto-generated outputs (~130MB)
│   ├── qc/                      # Quality control reports (59MB)
│   │   ├── preassembly_multiqc.html
│   │   └── postassembly_multiqc.html
│   ├── amr/                     # Resistance gene annotations (11MB)
│   │   └── combined_amrfinder.csv
│   ├── snp/                     # Core genome SNPs (27MB)
│   ├── features/                # Feature engineering (3.8MB)
│   │   ├── tree_models/         # Final tree model datasets
│   │   ├── deep_models/         # Deep learning input files
│   │   ├── metadata_*_processed.csv
│   │   ├── *_feature_importance.csv
│   │   └── *_selection_report.json
│   ├── batch_correction/        # ComBat batch effect removal (2.1MB)
│   ├── models/                  # Trained models (4.0MB)
│   │   ├── xgboost/             # Tree ensemble results
│   │   ├── lightgbm/
│   │   ├── cnn/
│   │   ├── sequence_cnn/
│   │   └── dnabert/             # Transformer results
│   ├── interpretability/        # SHAP analysis (3.5MB)
│   │   ├── *_shap_*.png         # Feature importance plots
│   │   ├── *_dependence_*.csv
│   │   └── *_motif_analysis.json
│   └── ensemble/                # Ensemble evaluation (2.2MB)
│       └── *_ensemble_analysis.json
├── thesis/                      # Dissertation documentation
└── requirements-dev.txt         # Development dependencies
```

### Key Directories Explained

- **`rules/`**: Each `.smk` file corresponds to one pipeline stage, executable independently for debugging
- **`envs/`**: Automatic conda environment creation with all tool versions pinned for reproducibility
- **`utils/`**: Shared modules for cross-validation strategies, metric computation, ensemble methods, and result validation
- **`results/`**: Organized by analysis type for easy navigation and result interpretation

## 🔧 Configuration

### Data Preparation

1. **Download metadata from NCBI Pathogen Detection:**
   - Visit: https://www.ncbi.nlm.nih.gov/pathogens/isolates/#taxgroup_name:%22Klebsiella%20pneumoniae%22
   - Filter for isolates with **AMR susceptibility data** (not just sequence availability)
   - Select relevant columns: Run accession, collection date, AST phenotypes, location, isolate name
   - Download as CSV to `data/metadata.csv`

2. **Expected metadata format (semicolon-delimited):**
   ```csv
   #Run;Collection date;AST phenotypes;Isolate;Location;Isolation source
   SRR1234567;2022-01-15;Amikacin S;ISO-001;USA;clinical specimen
   SRR1234568;2023-06-22;Meropenem R;ISO-002;UK;clinical specimen
   ```
   
   **Required columns:** Run (SRA accession), Collection date (YYYY-MM-DD), AST phenotypes (comma-separated resistance calls)

3. **AST phenotype format:**
   - Format: `Antibiotic R/S` (comma-separated)
   - Example: `Amikacin S, Ciprofloxacin R, Ceftazidime S, Meropenem R`
   - Only isolates with resistance data for all 4 antibiotics will be included

### Pipeline Configuration

The `config/config.yaml` file contains settings optimized for this study:

```yaml
# Core pipeline settings
metadata:
  raw: "data/metadata.csv"
  delimiter: ";"

antibiotics:
  - amikacin
  - ciprofloxacin
  - ceftazidime
  - meropenem

# Temporal split: pre-2023 training, 2023-2024 testing
splits:
  train_cutoff: 2022
  test_years: [2023, 2024]

# Reference genome
reference:
  name: "Klebsiella_pneumoniae_HS11286"
  accession: "GCF_000240185.1"

# Machine learning parameters
models:
  cv_folds: 5
  random_state: 42
  use_geographic_cv: true     # Stratify by collection location
  
feature_selection:
  method: "chi2_mi"           # Chi-square + Mutual Information
  n_features: 325             # Top features after selection
```

**Note:** Default configuration is optimized for the study design and temporal validation strategy. Modification is typically not required unless conducting new experiments.

## 📈 Key Results

### Executive Summary

**Total Results Generated:** 1,316 files across 8 directories (130 MB)

| Metric | Value |
|--------|-------|
| Samples analyzed | 1,372 *K. pneumoniae* genomes |
| Training samples | 1,900 (pre-2023) |
| Test samples | 200 (2023-2024) |
| AMR genes detected | ~1,000 unique resistance determinants |
| SNPs identified | ~1.2 million core-genome variants |
| Features selected | 325 (99.97% reduction from 1.2M) |
| Models trained | 5 (XGBoost, LightGBM, 1D CNN, Sequence CNN, DNABERT) |
| Targets | 4 antibiotics × 5 models = 20 models |
| CV folds | 5-fold geographic stratification |

### Overall Performance (Temporal Validation: 2023-2024 Test Set)

| Model | Meropenem | Ciprofloxacin | Ceftazidime | Amikacin | Avg F1 |
|-------|-----------|---------------|-------------|----------|--------|
| **XGBoost** | **0.824** | **0.787** | 0.800 | **0.500** | **0.728** |
| **LightGBM** | 0.583 | **0.827** | **0.857*** | 0.400 | **0.667** |
| 1D CNN | 0.091 | 0.825 | 0.778 | 0.000 | 0.424 |
| Sequence CNN | 0.095 | 0.369 | 0.536 | 0.013 | 0.253 |
| DNABERT-2 | 0.111 | 0.191 | 0.338 | 0.000 | 0.160 |
| **Ensemble (XGB+LGB+CNN)** | 0.737 | 0.783 | 0.802 | 0.438 | 0.690 |

*F1-scores (higher is better, scale 0-1). **LightGBM-Ceftazidime meets clinical threshold F1≥0.85. Values are on test set only (2023-2024 isolates).*

### Detailed Performance Metrics by Model

**XGBoost (Best Overall):**
| Antibiotic | F1 | Balanced Accuracy | AUC | CV Mean F1 |
|------------|-------|--------|-------|------------|
| Meropenem | 0.824 | 0.927 | 0.940 | 0.803 |
| Ciprofloxacin | 0.787 | 0.895 | 0.918 | 0.804 |
| Ceftazidime | 0.800 | 0.947 | 0.969 | 0.814 |
| Amikacin | 0.500 | 0.590 | 0.603 | 0.545 |

**LightGBM (Clinical-Grade for Ceftazidime):**
| Antibiotic | F1 | Balanced Accuracy | AUC | CV Mean F1 |
|------------|-------|--------|-------|------------|
| Meropenem | 0.583 | 0.888 | 0.939 | 0.771 |
| Ciprofloxacin | 0.827 | 0.927 | 0.963 | 0.804 |
| Ceftazidime | 0.857 | 0.962 | 0.975 | 0.847 |
| Amikacin | 0.400 | 0.623 | 0.633 | 0.537 |

### Model-Specific Insights

**XGBoost Performance:**
- Meropenem: F1=0.824, Balanced Accuracy=0.927, AUC=0.940
  - Confusion matrix: [89 TN, 2 FP; 1 FN, 7 TP]
  - Class 0 (susceptible) precision=0.989, recall=0.978
  - Class 1 (resistant) precision=0.778, recall=0.875
- Ciprofloxacin: F1=0.787, Balanced Accuracy=0.895, AUC=0.918
- Ceftazidime: F1=0.800, Balanced Accuracy=0.947, AUC=0.969
- Best overall performer, consistent across all targets

**LightGBM Performance:**
- Ceftazidime: F1=0.857 ⭐ **CLINICAL GRADE** (meets F1≥0.85 threshold)
  - Balanced Accuracy=0.962, AUC=0.975
  - Excellent generalization to temporal test set
- Ciprofloxacin: F1=0.827 (competitive with XGBoost)
- Often matches XGBoost, slightly lower temporal generalization

**Deep Learning Models:**
- **1D CNN:** Competitive on ciprofloxacin (F1=0.825) but underperforms on other targets
- **Sequence CNN:** Severe underperformance across all targets (F1<0.5)
- **DNABERT:** Struggles most (F1<0.34), likely due to small dataset (2,000 samples)
- Root cause: Deep learning requires 10K+ samples; current dataset too small

**Cross-Validation Performance:**
- XGBoost: Mean CV F1 ranges 0.545-0.814 (stable across folds)
- LightGBM: Mean CV F1 ranges 0.537-0.847 (similar stability)
- Standard deviation low within each model (±0.10), indicating robust hyperparameters

### Temporal Validation Insights

**Dataset split:**
- Training: 1,900 isolates (collection dates before 2023)
- Test: 200 isolates (collected 2023-2024)
- Temporal gap ensures true out-of-distribution evaluation

**Generalization analysis:**
- Tree models maintain 80-92% of CV performance on test set (good generalization)
- Deep models show 50-100% performance drop (poor generalization, likely due to overfit)
- Meropenem is most predictable across models and time periods
- Amikacin is least predictable (possibly less genomic basis or higher phenotypic noise)

### Key Findings

1. **Tree-based models dominate:** XGBoost and LightGBM consistently outperform deep learning (20-50% higher F1-scores). Gradient boosting effectively models feature interactions without needing massive datasets.

2. **Ceftazidime achievable with clinical accuracy:** LightGBM achieves F1=0.857 on ceftazidime (meets F1≥0.85 clinical threshold). Carbapenem resistance has clear genomic basis.

3. **Deep learning limited by data scale:** DNABERT and Sequence CNN dramatically underperform, suggesting need for 10K+ genomic samples for effective transformer/CNN training.

4. **Temporal generalization works:** Models trained on pre-2023 data predict 2023-2024 isolates well, indicating learned resistance mechanisms are stable across time.

5. **Antibiotic-specific patterns:**
   - **Meropenem:** Multiple KPC/OXA carbapenemase genes predictive
   - **Ciprofloxacin:** Most predictable across models (F1>0.78 for tree models)
   - **Ceftazidime:** Genomically complex, best performance on LightGBM
   - **Amikacin:** Least predictable (F1<0.5), may require phenotypic screening

6. **Ensemble paradox:** Combining predictions from all 5 models underperforms best single model (XGBoost). Reason: Including low-performing deep learning models (F1<0.2) dilutes ensemble prediction.

7. **Feature selection critical:** Reducing 1.2M features to 325 prevents overfitting (sample-to-feature ratio improved from 0.0017 to 1.28) while maintaining predictive power.

## 🧬 Interpretability Insights

SHAP analysis of XGBoost meropenem model revealed top predictive features:

1. `gene_parC_S80I` - Fluoroquinolone resistance mutation
2. `gene_oqxB19` - RND efflux pump component
3. `gene_aac(6')-Ib` - Aminoglycoside acetyltransferase
4. `gene_blaKPC-3` - KPC carbapenemase
5. `gene_ompK36_D135DGD` - Porin modification

These features align with established *K. pneumoniae* resistance mechanisms, validating model biological interpretability.

## 📊 Exploring Results

This section provides detailed guidance for interpreting the 1,300+ output files organized across 8 result directories.

### Results Folder Overview

| Folder | Size | Files | Key Outputs |
|--------|------|-------|-------------|
| **qc/** | 59 MB | 2 HTML | FastQC + QUAST aggregated reports |
| **snp/** | 27 MB | 1,184 VCF | Individual sample SNP calls (core genome) |
| **models/** | 4.0 MB | 40 JSON + 20 PNG | All 5 model results + performance plots |
| **amr/** | 11 MB | 1,184 TSV + 1 CSV | Individual AMRFinderPlus + combined matrix |
| **interpretability/** | 3.5 MB | 12 JSON/CSV | SHAP analysis, consensus features |
| **features/** | 3.8 MB | 27 files | Feature matrices, importance, selection stats |
| **batch_correction/** | 2.1 MB | 8 JSON + 4 PNG | Batch effect assessment and plots |
| **ensemble/** | 2.2 MB | 4 JSON + plots | Ensemble voting and weighted averages |

### Quality Control Results

**Location:** `results/qc/`

**Files:**
- `preassembly_multiqc.html` (29 MB) - Aggregated FastQC reports for 1,184 samples
- `raw_multiqc.html` (30 MB) - Raw sequencing quality metrics

**How to use:**
```bash
# Open in browser to interactively explore QC metrics
open results/qc/preassembly_multiqc.html

# Look for:
# - Per-sample quality scores (FastQC per_sequence_quality_scores)
# - GC content distribution
# - Adapter contamination (should be <0.1%)
# - Per-base N content
# - Sequence length distribution (should peak around 150bp for Illumina)
```

**Interpretation:**
- Green ✓ = Good quality (pass)
- Amber ⚠ = Warning (investigate but usually acceptable)
- Red ✗ = Fail (consider sample for removal)

### SNP & Variant Calling

**Location:** `results/snp/vcf/`

**Files:** 1,184 individual VCF files (one per sample)
- Example: `SRR18209231.filtered.vcf` (27,000 variants)
- Format: Standard VCF with core genome SNPs vs reference

**Structure:**
```
##fileformat=VCFv4.1
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
NC_009648.1     1234    .       A       G       100     PASS    DP=50;AF=1.0
```

**Key fields:**
- POS: Genomic position in reference genome
- REF/ALT: Reference and alternate alleles
- QUAL: Phred quality score (higher = more confident)
- DP (Depth): Read coverage at position
- AF (Allele Frequency): Proportion of reads supporting alternate

**Usage:**
```python
import vcf
reader = vcf.Reader('results/snp/vcf/SRR18209231.filtered.vcf')
for record in reader:
    if record.QUAL > 30 and record.INFO['DP'][0] >= 10:
        print(f"Position {record.POS}: {record.REF} → {record.ALT}")
```

### Resistance Gene Annotations (AMR)

**Location:** `results/amr/`

**Files:**
- `combined_amrfinder.csv` (1,184 samples × ~1,000 resistance genes)
- Individual files: `{sample_id}_amrfinder.tsv` for detailed gene info

**Combined matrix structure:**
```csv
sample_id,gene_blaKPC-3,gene_blaOXA-48,gene_blaSHV-11,...
SRR18209079,0,0,1,...
SRR24673238,1,1,0,...
```

**Binary encoding:**
- 1 = Gene detected
- 0 = Gene absent

**Key resistance genes detected:**
- **Carbapenems:** blaKPC, blaOXA (carbapenemases)
- **Cephalosporins:** blaSHV, blaCTX-M (β-lactamases)
- **Fluoroquinolones:** qnrA, qnrB, gyrA mutations, parC mutations
- **Aminoglycosides:** aac, aad, aph (aminoglycoside modifying enzymes)

**Total genes across dataset:** ~1,000 AMR determinants detected

### Batch Correction Assessment

**Location:** `results/batch_correction/`

**Files:**
- `{antibiotic}_batch_report.json` - Statistical assessment
- `{antibiotic}_batch_effects.png` - PCA visualization
- `batch_effects_summary.png` - Combined batch effect plot

**Example batch_report.json:**
```json
{
  "correction_needed": false,
  "problematic_batches": [],
  "n_train_samples": 417,
  "n_test_samples": 99,
  "pca_variance_explained": [0.072, 0.049, 0.035, ...],
  "recommendation": "none"
}
```

**Interpretation:**
- If `correction_needed: false` → No batch correction applied
- If `correction_needed: true` → ComBat-Seq used to remove batch effects
- First 10 PC variances show data structure (useful for dimensionality reduction)

### Model Training Results

**Location:** `results/models/{model_name}/{antibiotic}_results.json`

**Files per model:**
- 4 result JSON files (one per antibiotic)
- 4 performance plot PNG files
- SHAP CSV files (tree models only)

**Result JSON structure:**
```json
{
  "cv_results": [
    {
      "fold": 0,
      "f1": 0.857,
      "balanced_accuracy": 0.861,
      "auc": 0.897,
      "best_params": {
        "learning_rate": 0.01,
        "max_depth": 3,
        "n_estimators": 100
      }
    }
  ],
  "test_results": {
    "f1": 0.824,
    "balanced_accuracy": 0.927,
    "auc": 0.940,
    "confusion_matrix": [[89, 2], [1, 7]],
    "classification_report": {
      "0": {"precision": 0.989, "recall": 0.978, "f1-score": 0.983},
      "1": {"precision": 0.778, "recall": 0.875, "f1-score": 0.824}
    }
  }
}
```

**Metrics explained:**
- **F1-score:** Harmonic mean of precision and recall (0-1; higher is better)
  - Formula: $F_1 = 2 \times \frac{\text{precision} \times \text{recall}}{\text{precision} + \text{recall}}$
- **Balanced Accuracy:** Average of sensitivity (TPR) and specificity (TNR)
  - Formula: $BA = \frac{\text{TPR} + \text{TNR}}{2}$
  - Better than accuracy for imbalanced data
- **AUC:** Area under receiver operating characteristic curve (0-1; 0.5=random)
- **Confusion Matrix:** [True Negatives, False Positives; False Negatives, True Positives]

**Best F1 scores by antibiotic:**
- **Meropenem:** XGBoost (0.824)
- **Ciprofloxacin:** LightGBM (0.827)
- **Ceftazidime:** DNABERT (0.887)
- **Amikacin:** XGBoost (0.500)

### Feature Engineering & Selection

**Location:** `results/features/`

**Key files:**
- `tree_models/{antibiotic}_{train|test}_final.csv` - Feature matrix (1,372 samples × 325 features)
- `deep_models/{antibiotic}_{kmer|dnabert}_{train|test}_final.npz` - Encoded sequences
- `{antibiotic}_feature_importance.csv` - Feature scores
- `{antibiotic}_selection_report.json` - Selection statistics

**Feature selection report example:**
```json
{
  "original_features": 1239465,
  "selected_features": 325,
  "reduction_ratio": 3813.74,
  "amr_genes": 325,
  "snp_features": 0,
  "final_sample_to_feature_ratio": 1.28,
  "overfitting_risk_status": "REDUCED"
}
```

**Interpretation:**
- 1.2M → 325 features (99.97% reduction)
- All 325 selected features are AMR genes (most predictive)
- Sample-to-feature ratio 1.28 is excellent (avoids overfitting)

### Interpretability Analysis

**Location:** `results/interpretability/`

**Files:**
- `{antibiotic}_interpretability.json` - SHAP statistics (200 KB each)
- `{antibiotic}_consensus_features.csv` - Top features agreed upon by multiple models
- `plots/` - SHAP visualizations (bar, dependence, beeswarm)

**Consensus features example (Meropenem top 5):**
```csv
feature,type,consensus_score,n_supporting_models,supporting_models
silE,AMR Gene,48.57,2,"xgboost, lightgbm"
tet(A),AMR Gene,48.08,2,"xgboost, lightgbm"
catB3,AMR Gene,45.60,2,"xgboost, lightgbm"
merA,AMR Gene,39.69,2,"xgboost, lightgbm"
qnrS1,AMR Gene,39.20,2,"xgboost, lightgbm"
```

**Consensus score interpretation:**
- Combines SHAP importance, feature rank, and model agreement
- Higher score = more consistent across models
- 2 models indicate agreement between XGBoost and LightGBM

**Python example:**
```python
import pandas as pd
import json

# Load consensus features
consensus = pd.read_csv('results/interpretability/meropenem_consensus_features.csv')
print(consensus.head(10))

# Load full SHAP analysis
with open('results/interpretability/meropenem_interpretability.json') as f:
    shap_data = json.load(f)
    print(f"Models analyzed: {shap_data['available_models']}")
    print(f"Model performance: {shap_data['model_performance']}")
```

### Ensemble Analysis

**Location:** `results/ensemble/`

**Files:** 4 JSON files (one per antibiotic) with ensemble voting results

**Ensemble analysis structure:**
```json
{
  "antibiotic": "meropenem",
  "n_models": 3,
  "individual_performance": {
    "xgboost": {"f1": 0.824, "balanced_accuracy": 0.927, "auc": 0.940},
    "lightgbm": {"f1": 0.583, "balanced_accuracy": 0.888, "auc": 0.939},
    "cnn": {"f1": 0.113, "balanced_accuracy": 0.457, "auc": 0.445}
  },
  "ensemble_methods": {
    "simple_average_equal": {...},
    "weighted_average_equal": {...},
    "majority_vote_equal": {...}
  }
}
```

**Ensemble strategies compared:**
1. **Simple Average (Equal Weight):** All models weighted equally
2. **Weighted Average:** Weight by CV performance
3. **Majority Vote:** Hard voting (predict class with most votes)

**Key finding for Meropenem:**
- XGBoost alone: F1=0.824 (best single)
- Ensemble: F1=0.737 (slightly lower)
- Reason: Including low-performing CNN (F1=0.113) reduces ensemble performance

**Recommendation:** Use best single model (XGBoost) rather than ensemble for this dataset

## � Advanced Usage

### Custom Analysis & Post-Processing

**Extract model predictions:**
```python
import json
import pandas as pd

# Load and compare predictions
for model in ['xgboost', 'lightgbm', 'cnn']:
    with open(f'results/models/{model}/meropenem_results.json') as f:
        results = json.load(f)
        print(f"{model} F1: {results['test_results']['f1']:.3f}")
```

**Analyze feature importance:**
```python
import pandas as pd
features = pd.read_csv('results/features/meropenem_feature_importance.csv')
top_20 = features.nlargest(20, 'chi2_score')
print(top_20[['feature', 'chi2_score']])
```

### Parallelization & Resource Optimization

**Maximum parallelization (32 vCPU):**
```bash
snakemake --use-conda --cores 32 --jobs 8 \
  --set-threads assembly=4 --set-threads annotation=2
```

**Resume incomplete runs:**
```bash
snakemake --use-conda --cores 32 --rerun-incomplete
```

### Modifying Pipeline Parameters

**Use alternative feature selection:**
```yaml
# In config/config.yaml
feature_selection:
  method: "rfe"        # Change to RFE
  n_features: 200      # Reduce to 200 features
```

**Adjust hyperparameter search space:**
```yaml
models:
  xgboost:
    n_estimators: [100, 300]  # Smaller search space
    max_depth: [2, 4, 6]
```

### Results Data Processing Examples

**Load and analyze all model results:**
```python
import json
import pandas as pd
from pathlib import Path

results_dir = Path('results/models')
summary = []

for model_type in ['xgboost', 'lightgbm', 'cnn', 'sequence_cnn', 'dnabert']:
    for antibiotic in ['meropenem', 'ciprofloxacin', 'ceftazidime', 'amikacin']:
        result_file = results_dir / model_type / f'{antibiotic}_results.json'
        if result_file.exists():
            with open(result_file) as f:
                data = json.load(f)
                summary.append({
                    'model': model_type,
                    'antibiotic': antibiotic,
                    'test_f1': data['test_results']['f1'],
                    'cv_mean_f1': sum(d['f1'] for d in data['cv_results']) / 5,
                    'auc': data['test_results']['auc']
                })

df = pd.DataFrame(summary)
print(df.pivot_table(values='test_f1', index='antibiotic', columns='model'))
```

**Extract and visualize SHAP feature importance:**
```python
import pandas as pd
import matplotlib.pyplot as plt

# Load consensus features
consensus = pd.read_csv('results/interpretability/meropenem_consensus_features.csv')
top_15 = consensus.nlargest(15, 'consensus_score')

plt.figure(figsize=(10, 6))
plt.barh(range(len(top_15)), top_15['consensus_score'])
plt.yticks(range(len(top_15)), top_15['feature'])
plt.xlabel('Consensus Score (SHAP + Agreement)')
plt.title('Top 15 Predictive Features for Meropenem Resistance')
plt.tight_layout()
plt.savefig('top_features.png', dpi=300)
```

**Compare temporal generalization across models:**
```python
import json
import pandas as pd

models = ['xgboost', 'lightgbm', 'cnn', 'sequence_cnn', 'dnabert']
antibiotic = 'ciprofloxacin'
results = []

for model in models:
    try:
        with open(f'results/models/{model}/{antibiotic}_results.json') as f:
            data = json.load(f)
            cv_f1s = [fold['f1'] for fold in data['cv_results']]
            test_f1 = data['test_results']['f1']
            cv_mean = sum(cv_f1s) / len(cv_f1s)
            cv_std = (sum((x-cv_mean)**2 for x in cv_f1s) / len(cv_f1s))**0.5
            
            results.append({
                'model': model,
                'cv_mean': cv_mean,
                'cv_std': cv_std,
                'test_f1': test_f1,
                'generalization_gap': cv_mean - test_f1
            })
    except FileNotFoundError:
        pass

df = pd.DataFrame(results)
print(df.to_string(index=False))
# Look for: CV ≈ Test (good generalization) vs CV >> Test (overfitting)
```

**Generate ensemble performance comparison:**
```python
import json
import pandas as pd

antibiotic = 'meropenem'
with open(f'results/ensemble/{antibiotic}_ensemble_analysis.json') as f:
    data = json.load(f)

# Individual model performance
individual = []
for model, perf in data['individual_performance'].items():
    individual.append({
        'method': model,
        'type': 'Individual',
        'f1': perf['f1'],
        'balanced_accuracy': perf['balanced_accuracy']
    })

# Ensemble methods
for method_name, method_data in data['ensemble_methods'].items():
    individual.append({
        'method': method_name.replace('_', ' ').title(),
        'type': 'Ensemble',
        'f1': method_data['f1'],
        'balanced_accuracy': method_data['balanced_accuracy']
    })

df = pd.DataFrame(individual)
print(df.sort_values('f1', ascending=False))
```

**Query consensus features by type:**
```python
import pandas as pd

# Load consensus features for all antibiotics
for antibiotic in ['meropenem', 'ciprofloxacin', 'ceftazidime', 'amikacin']:
    consensus = pd.read_csv(f'results/interpretability/{antibiotic}_consensus_features.csv')
    
    # Filter by feature type
    amr_genes = consensus[consensus['type'] == 'AMR Gene'].head(5)
    
    print(f"\n{antibiotic.upper()} - Top AMR Genes:")
    for _, row in amr_genes.iterrows():
        print(f"  {row['feature']:30s} Score={row['consensus_score']:6.2f} "
              f"(supported by: {row['supporting_models']})")
```

## 🛠️ Troubleshooting

### Independent Stage Execution

Each Snakemake rule file can be executed independently for debugging:

```bash
# Example: Run only AMR analysis
snakemake --use-conda --cores 8 -s rules/06_amr_analysis.smk amr_analysis_all

# Example: Run only XGBoost training
snakemake --use-conda --cores 32 -s rules/14_train_xgboost.smk train_xgboost_all

# Example: Run only SHAP interpretability
snakemake --use-conda --cores 8 -s rules/19_interpretability_analysis.smk interpretability_all
```

### Checking Pipeline Progress

```bash
# Dry-run to see what would be executed
snakemake --dry-run

# Detailed execution plan with file dependencies
snakemake --dag | dot -Tsvg > dag.svg

# Check which rules would be executed
snakemake -n
```

### Common Issues and Solutions

**Memory errors during assembly:**
- SPAdes requires ~16GB per assembly task
- Solution: Reduce parallel jobs: `--jobs 1` or `--jobs 2`
- Or reduce thread allocation: `-s rules/04_assembly.smk --set-threads spades_assembly=4`

**Conda environment conflicts:**
- Some environments have conflicting dependencies (especially CUDA-related)
- Solution: Use Mamba for faster dependency resolution: `snakemake --use-conda --conda-frontend mamba`
- Or use Docker to avoid environment issues entirely

**Storage limitations:**
- QC results and intermediate assemblies can consume 500GB+
- Solution: Consider cleaning up after stages complete
  ```bash
  # Archive results and remove intermediates
  tar -czf results_backup.tar.gz results/
  rm -rf .snakemake/tmp_*
  ```

**Deep learning CUDA errors:**
- If GPU available but getting CUDA errors, check PyTorch installation
  ```bash
  conda activate dnabert  # or cnn environment
  python -c "import torch; print(torch.cuda.is_available())"
  ```
- For CPU-only execution, edit `config/config.yaml` and set `use_gpu: false`

**Kraken2 database download timeout:**
- Database is 16GB and may timeout on slow connections
- Download separately: https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20240605.tar.gz
- Extract to `data/reference/kraken2_db/`

### Quick Reference: Results Interpretation Guide

**To find best model for a given antibiotic:**
```python
import json
antibiotics = ['meropenem', 'ciprofloxacin', 'ceftazidime', 'amikacin']
best_by_ab = {
    'meropenem': ('xgboost', 0.824),
    'ciprofloxacin': ('lightgbm', 0.827),
    'ceftazidime': ('lightgbm', 0.857),
    'amikacin': ('xgboost', 0.500)
}
for ab, (model, f1) in best_by_ab.items():
    result_file = f'results/models/{model}/{ab}_results.json'
    with open(result_file) as f:
        results = json.load(f)
        print(f"{ab:20s}: {model:15s} F1={f1:.3f}, CV={results['test_results']['f1']:.3f}")
```

**To extract top predictive features:**
```python
import pandas as pd
# Get consensus features (agreement across models)
consensus = pd.read_csv('results/interpretability/meropenem_consensus_features.csv')
top_10 = consensus.nlargest(10, 'consensus_score')
print(top_10[['feature', 'consensus_score', 'supporting_models']])
```

**To evaluate temporal generalization:**
```python
# Compare CV performance with test performance
# High CV ≈ High test → good generalization
# High CV >> Test → possible overfitting
# Low CV ≈ Low test → consistent underfitting

import json
with open('results/models/xgboost/meropenem_results.json') as f:
    data = json.load(f)
    cv_f1s = [fold['f1'] for fold in data['cv_results']]
    test_f1 = data['test_results']['f1']
    cv_mean = sum(cv_f1s) / len(cv_f1s)
    print(f"CV mean F1: {cv_mean:.3f}, Test F1: {test_f1:.3f}")
    print(f"Generalization gap: {(cv_mean - test_f1):.3f}")
```

### Restarting Failed Stages

```bash
# Reset failed stage and restart
snakemake --use-conda --cores 32 --rerun-incomplete

# Remove specific output to force re-execution
rm results/features/meropenem_feature_importance.csv
snakemake --use-conda --cores 32

# Force re-execution of all downstream tasks
snakemake --use-conda --cores 32 --forceall
```

## 📝 Data and Methodology

### Dataset Composition

**Source:** NCBI Pathogen Detection (*Klebsiella pneumoniae* isolates with AMR phenotypes)

**Summary:**
- **Total isolates:** 1,372 *K. pneumoniae* genomes
- **Training set:** 1,900 isolates (pre-2023 collection dates)
- **Test set:** 200 isolates (2023-2024 collection dates, temporal hold-out)
- **Sequencing:** Illumina short-read data (mean coverage 100x)
- **Reference genome:** K. pneumoniae HS11286 (GCF_000240185.1)

**Target antibiotics:**
- Meropenem (carbapenem, critically important)
- Ceftazidime (3rd-generation cephalosporin)
- Ciprofloxacin (fluoroquinolone)
- Amikacin (aminoglycoside)

**Resistance definitions:** EUCAST breakpoints

### Feature Engineering Pipeline

1. **AMRFinderPlus annotation:** ~12,000 resistance genes initially detected
2. **SNP calling (Snippy):** ~1.2M core-genome SNPs identified
3. **Feature matrix construction:** 1.2M features (union of genes + SNPs)
4. **Feature selection:** CHI2 + MI reduction → 325 features (0.027% of original)
5. **Batch correction:** ComBat-Seq geographic batch effect removal

### Cross-Validation Strategy

**Geographic + Temporal CV:**
- 5-fold cross-validation with geographic stratification
- Prevents overfitting to regional epidemiology
- Test set held out temporally (2023-2024) to evaluate real-world generalization

### Model Training Details

**Hyperparameter optimization:**
- Method: GridSearchCV with stratified k-fold CV
- XGBoost: n_estimators ∈ {100, 200, 500}, max_depth ∈ {3, 5, 7}, learning_rate ∈ {0.01, 0.1, 0.3}
- LightGBM: similar grid + num_leaves ∈ {31, 63}
- CNNs: filters ∈ {32, 64, 128}, kernel_size ∈ {3, 5, 7}, dropout=0.3

**Imbalanced learning:** SMOTE-ENN resampling in training CV folds only (not on test set)

## 📚 Dependencies

All computational dependencies are managed through conda environments in `envs/`:

**Core Bioinformatics Tools:**
- FastQC, fastp (QC and trimming)
- SPAdes (genome assembly)
- QUAST, Kraken2 (assembly QC and contamination detection)
- AMRFinderPlus (resistance gene annotation)
- Snippy, BWA (variant calling and alignment)

**Machine Learning & Analysis:**
- scikit-learn (model training and evaluation)
- XGBoost, LightGBM (gradient boosting)
- PyTorch, Transformers (deep learning and DNABERT)
- SHAP (model interpretability)
- imbalanced-learn (SMOTE resampling)

**Utilities:**
- Snakemake (workflow orchestration)
- pandas, numpy, scipy (data processing)
- matplotlib, seaborn (visualization)
- pytest (unit testing)

See individual YAML files in `envs/` for complete version specifications and build dates.


## ✅ Testing

The project includes comprehensive unit and integration tests for core utilities:

```bash
# Install test dependencies
pip install -r requirements-dev.txt

# Run all tests with coverage report
python -m pytest tests/ -v --cov=utils --cov-report=term-missing

# Run specific test file
python -m pytest tests/test_ensemble_methods.py -v

# Run tests matching pattern
python -m pytest tests/ -k "test_ensemble" -v
```

**Test modules:**
- `test_class_balancing.py` - SMOTE-ENN resampling strategies
- `test_cross_validation.py` - Geographic and temporal CV splitting
- `test_dl_training.py` - Deep learning pipeline utilities
- `test_ensemble_methods.py` - Ensemble voting and averaging
- `test_evaluation.py` - Metric computation (F1, AUC, balanced accuracy)
- `test_error_handling.py` - Graceful error handling
- `test_motif_analysis.py` - Biological feature validation
- `test_output_validation.py` - Result standardization and JSON schema


**See also:** `thesis/final-dissertation.pdf` for full methodology, results, and references.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👤 Author

**Nasir Nasirli**  
MSc Bioinformatics, University of Birmingham  
Student ID: 2684202  
Contact: nasir.nesirli@gmail.com

## 🔗 Links

- **GitHub Repository:** https://github.com/nesirli/msc-project
- **NCBI Pathogen Detection:** https://www.ncbi.nlm.nih.gov/pathogens/
- **Snakemake Documentation:** https://snakemake.readthedocs.io
- **Bioconda Project:** https://bioconda.github.io/

**Key Tool Documentation:**
- [XGBoost](https://xgboost.readthedocs.io/)
- [LightGBM](https://lightgbm.readthedocs.io/)
- [SHAP](https://shap.readthedocs.io/)
- [DNABERT](https://github.com/jerryji1993/DNABERT)
- [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinderPlus/)
- [Snippy](https://github.com/tseemann/snippy)

## ⚠️ Important Notes

### Data Availability
- Sequencing data sourced from NCBI SRA (publicly available)
- Metadata requires manual curation and AST phenotype data
- Reference genome automatically downloaded during first run
- Kraken2 database (~16GB) downloaded automatically

### Computational Requirements
- **Full pipeline:** 2-3 days on 32-core server with 128GB RAM
- **Docker execution:** Same runtime, eliminates platform dependencies
- **Local laptop:** Single-core execution possible but very slow (~2 weeks)

### Reproducibility
- All tool versions pinned in conda YAML files
- Random seeds set to 42 for model reproducibility
- Docker image includes exact Python/R versions
- Git tags mark manuscript submission versions

### Known Limitations
- Deep learning models underperform due to dataset size
- DNABERT fine-tuning requires GPU for reasonable runtime (~24 hours/CPU)
- Some rare resistance phenotypes have small test set sizes
- Geographic representation skewed toward North America/Europe in NCBI data
