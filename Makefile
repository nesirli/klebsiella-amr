# Klebsiella pneumoniae AMR prediction pipeline
# Orchestrated with GNU Make + conda.
#
# macOS ships BSD make by default; this Makefile uses GNU Make features
# (group targets, pattern rules). On macOS install gmake and run 'gmake'.

ifeq ($(findstring .,$(MAKE_VERSION)),)
  $(error This Makefile requires GNU Make. On macOS: brew install make && gmake)
endif

.DELETE_ON_ERROR:

# Auto-generated includes ------------------------------------------------------
# Only include generated config if the amr environment exists; otherwise
# 'make setup' would try to build config.mk before the environment exists.
CONDA_BASE := $(shell conda info --base 2>/dev/null)
AMR_ENV := $(wildcard $(CONDA_BASE)/envs/amr)
ifneq ($(AMR_ENV),)
  -include results/metadata/config.mk
  ifneq ($(wildcard results/metadata/samples.mk),)
    include results/metadata/samples.mk
  endif
endif

# Conda runners ----------------------------------------------------------------
RUN_AMR     := conda run --no-capture-output -n amr
RUN_BIOINFO := conda run --no-capture-output -n bioinfo

# Settings ---------------------------------------------------------------------
MAX_SAMPLES  := -1
BATCH_SIZE   ?= 1
DOWNLOAD_JOBS ?= 4

TUNE_TRIALS_XGB := 30
TUNE_TRIALS_LGB := 30
TUNE_TRIALS_NN  := 20
TUNE_SPLITS     := 3

# Directories ------------------------------------------------------------------
RESULTS_DIR     := results
DATA_DIR        := data
READS_DIR       := $(DATA_DIR)/reads
TRIMMED_DIR     := $(DATA_DIR)/trimmed
QC_DIR          := $(DATA_DIR)/qc
KRAKEN_DIR      := $(DATA_DIR)/kraken2
DOWNSAMPLED_DIR := $(DATA_DIR)/downsampled
ASSEMBLY_DIR    := $(DATA_DIR)/assembly
QUAST_DIR       := $(DATA_DIR)/quast
AMR_DIR         := $(DATA_DIR)/amr
CLEAN_DIR       := $(DATA_DIR)/cleaned
MULTIQC_DIR     := $(RESULTS_DIR)/multiqc
FEATURES_DIR    := $(RESULTS_DIR)/features
SEQUENCES_DIR   := $(RESULTS_DIR)/sequences
MODELS_DIR      := $(RESULTS_DIR)/models
REPORT_DIR      := $(RESULTS_DIR)/report

# Pipeline modules -------------------------------------------------------------
include modules/config.mk
include modules/metadata.mk
include modules/pipeline.mk
include modules/features.mk
include modules/models.mk
include modules/batch.mk

# Directory creation -----------------------------------------------------------
$(READS_DIR) $(TRIMMED_DIR) $(QC_DIR) $(KRAKEN_DIR) $(QUAST_DIR) $(CLEAN_DIR):
	mkdir -p $@

$(DOWNSAMPLED_DIR) $(ASSEMBLY_DIR):
	mkdir -p $@

$(AMR_DIR):
	mkdir -p $@

$(FEATURES_DIR) $(SEQUENCES_DIR):
	mkdir -p $@

$(MODELS_DIR)/xgboost $(MODELS_DIR)/lightgbm $(MODELS_DIR)/nn $(MODELS_DIR)/dnabert:
	mkdir -p $@

$(MULTIQC_DIR) $(REPORT_DIR):
	mkdir -p $@

# Phony targets ----------------------------------------------------------------
.PHONY: all setup test metadata models dnabert multiqc report clean \
        tune process-samples _download-all _process-samples analyze

all:
	@if [ -z "$(SAMPLES)" ]; then \
		$(MAKE) --no-print-directory metadata MAX_SAMPLES=$(MAX_SAMPLES) && \
		$(MAKE) --no-print-directory all MAX_SAMPLES=$(MAX_SAMPLES); \
	else \
		$(MAKE) --no-print-directory report; \
	fi

setup:
	mamba env create -f envs/env-ml.yml -n amr
	mamba env create -f envs/env-bio.yml -n bioinfo

test:
	$(RUN_AMR) pytest tests/ -v

clean:
	rm -rf $(DATA_DIR) $(RESULTS_DIR) reference
