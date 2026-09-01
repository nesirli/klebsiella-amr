# Klebsiella pneumoniae AMR prediction pipeline
# Orchestrated with GNU Make + conda.
#
# macOS ships BSD make by default; this Makefile uses GNU Make features
# (group targets, pattern rules). On macOS install gmake and run 'gmake'.

ifeq ($(findstring .,$(MAKE_VERSION)),)
  $(error This Makefile requires GNU Make. On macOS: brew install make && gmake)
endif

.DELETE_ON_ERROR:

# Files reached through pattern-rule chains (fastp JSONs, kraken2 reports,
# AMRFinder TSVs, assemblies) count as intermediate, and make deletes
# intermediates on exit. .SECONDARY with no prerequisites keeps every one of
# them; .DELETE_ON_ERROR still removes the target of a failed recipe.
.SECONDARY:

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
# MAX_SAMPLES comes from config.yaml via the include above; '?=' leaves that
# value alone and only supplies a default when config.mk has not been built
# yet. A command-line MAX_SAMPLES=N overrides both.
MAX_SAMPLES    ?= -1
BATCH_SIZE     ?= 1

# XGBoost, LightGBM and torch all read OMP_NUM_THREADS and otherwise each
# claim every core, so under -j the model targets oversubscribe the machine.
# MODEL_THREADS comes from config.yaml; the default covers the first pass
# before config.mk has been generated.
export OMP_NUM_THREADS := $(or $(MODEL_THREADS),4)

# Lets any op the Apple Silicon GPU backend lacks fall back to the CPU rather
# than aborting a fit. Harmless on machines without MPS.
export PYTORCH_ENABLE_MPS_FALLBACK := 1

TUNE_TRIALS_XGB := 30
TUNE_TRIALS_LGB := 30
TUNE_TRIALS_NN  := 20
TUNE_SPLITS     := 3

# DNABERT-2 refits a 117M-parameter encoder for every fold of every trial, so
# a trial costs minutes rather than milliseconds. Its search is deliberately
# smaller than the others', and runs on a capped subsample (see
# tune_dnabert.py) to keep the whole thing to hours.
TUNE_TRIALS_DNABERT := 5
TUNE_SPLITS_DNABERT := 2
TUNE_MAX_TRAIN_DNABERT := 120

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

# Sample-list bootstrap --------------------------------------------------------
# The pipeline targets need SAMPLES, which is defined by the generated
# samples.mk included above. When it has not been built yet the variable is
# empty, so those targets build the metadata first and re-enter make to pick up
# the include. BOOTSTRAPPED marks the second pass: if SAMPLES is still empty
# there, the split genuinely produced nothing and we stop instead of recursing
# forever. Command-line variables propagate to sub-makes automatically.
#
# Usage: $(call need-samples,<real targets to build>)
define need-samples
@if [ -n "$(SAMPLES)" ]; then \
	$(MAKE) --no-print-directory $(1); \
elif [ -n "$(BOOTSTRAPPED)" ]; then \
	echo "error: no samples after 'make metadata' -- check that $(METADATA) has" >&2; \
	echo "       rows matching the train/test years in config.yaml" >&2; \
	exit 1; \
else \
	$(MAKE) --no-print-directory metadata && \
	$(MAKE) --no-print-directory $@ BOOTSTRAPPED=1; \
fi
endef

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
.PHONY: all setup test metadata models dnabert multiqc report _report clean force manifest app app-artifacts \
        tune process-samples _process-samples _process-one analyze

all:
	$(call need-samples,report)

setup:
	mamba env create -f envs/env-ml.yml -n amr
	mamba env create -f envs/env-bio.yml -n bioinfo

test:
	$(RUN_AMR) pytest tests/ -v

clean:
	rm -rf $(DATA_DIR) $(RESULTS_DIR) reference
