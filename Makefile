# Klebsiella pneumoniae AMR prediction pipeline
# Orchestrated with GNU Make + conda.

# macOS ships BSD make by default; this Makefile uses GNU Make features
# (group targets, pattern rules). On macOS install gmake and run 'gmake'.
ifeq ($(findstring .,$(MAKE_VERSION)),)
  $(error This Makefile requires GNU Make. On macOS: brew install make && gmake)
endif

# Delete partial outputs if a recipe fails, so incomplete files are not kept.
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

# Development sample cap. Override per run, e.g.:
#   make metadata MAX_SAMPLES=5
#   make metadata MAX_SAMPLES=20
#   make metadata MAX_SAMPLES=-1   # full dataset
MAX_SAMPLES := -1

# Optuna hyperparameter search settings. Tuning is run on the training split
# before final model training. Reduce trials for faster dev runs.
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

# Phony targets ----------------------------------------------------------------
# These point to stamp files so that the underlying data files can be treated as
# intermediates and removed automatically to save disk space.
.PHONY: all setup test metadata models dnabert multiqc report clean

all:
	@if [ -z "$(SAMPLES)" ]; then \
		$(MAKE) metadata MAX_SAMPLES=$(MAX_SAMPLES) && $(MAKE) all MAX_SAMPLES=$(MAX_SAMPLES); \
	else \
		$(MAKE) report; \
	fi

setup:
	mamba env create -f envs/env-ml.yml -n amr
	mamba env create -f envs/env-bio.yml -n bioinfo

test:
	$(RUN_AMR) pytest tests/ -v

# Config include ---------------------------------------------------------------
results/metadata/config.mk: config.yaml scripts/export_config.py
	mkdir -p results/metadata
	$(RUN_AMR) python3 scripts/export_config.py --config config.yaml --output $@

# Metadata ----------------------------------------------------------------------
metadata: results/metadata/.done

results/metadata/train.csv results/metadata/test.csv results/metadata/samples.mk results/metadata/samples.txt: results/metadata/.done
	@true

results/metadata/.done: config.yaml scripts/metadata.py results/metadata/config.mk
	mkdir -p results/metadata
	$(RUN_AMR) python3 scripts/metadata.py \
		--config config.yaml \
		--train-output results/metadata/train.csv \
		--test-output results/metadata/test.csv \
		--samples-output results/metadata/samples.txt \
		--max-samples $(MAX_SAMPLES)
	touch $@

# Reads ------------------------------------------------------------------------
$(READS_DIR)/.done: $(patsubst %,$(READS_DIR)/%_1.fastq.gz,$(SAMPLES)) | $(READS_DIR)
	@touch $@

$(READS_DIR)/%_1.fastq.gz $(READS_DIR)/%_2.fastq.gz &: scripts/download_reads.py | $(READS_DIR)
	$(RUN_AMR) python3 scripts/download_reads.py \
		--accession $* \
		--out1 $(READS_DIR)/$*_1.fastq.gz \
		--out2 $(READS_DIR)/$*_2.fastq.gz

$(READS_DIR):
	mkdir -p $@

# Trimming / fastp -------------------------------------------------------------
$(TRIMMED_DIR)/.done: $(patsubst %,$(QC_DIR)/%_fastp.json,$(SAMPLES)) | $(TRIMMED_DIR)
	@touch $@

$(TRIMMED_DIR)/%_1.fastq.gz $(TRIMMED_DIR)/%_2.fastq.gz $(QC_DIR)/%_fastp.json &: $(READS_DIR)/%_1.fastq.gz $(READS_DIR)/%_2.fastq.gz | $(TRIMMED_DIR) $(QC_DIR)
	$(RUN_BIOINFO) fastp \
		--in1 $(READS_DIR)/$*_1.fastq.gz \
		--in2 $(READS_DIR)/$*_2.fastq.gz \
		--out1 $(TRIMMED_DIR)/$*_1.fastq.gz \
		--out2 $(TRIMMED_DIR)/$*_2.fastq.gz \
		--json $(QC_DIR)/$*_fastp.json \
		--html $(QC_DIR)/$*_fastp.html \
		--qualified_quality_phred 20 \
		--length_required 50 \
		--thread $(FASTP_THREADS)

# Kraken2 ----------------------------------------------------------------------
$(KRAKEN_DB)/taxo.k2d:
	mkdir -p $(KRAKEN_DB)
	curl -L -o $(KRAKEN_DB)/k2_standard_08gb_20250402.tar.gz \
		https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz
	tar -xzf $(KRAKEN_DB)/k2_standard_08gb_20250402.tar.gz -C $(KRAKEN_DB)
	rm $(KRAKEN_DB)/k2_standard_08gb_20250402.tar.gz

$(KRAKEN_DIR)/.done: $(patsubst %,$(KRAKEN_DIR)/%_report.txt,$(SAMPLES)) | $(KRAKEN_DIR)
	@touch $@

$(KRAKEN_DIR)/%_report.txt: $(TRIMMED_DIR)/%_1.fastq.gz $(TRIMMED_DIR)/%_2.fastq.gz $(KRAKEN_DB)/taxo.k2d | $(KRAKEN_DIR)
	$(RUN_BIOINFO) kraken2 --db $(KRAKEN_DB) \
		--paired --gzip-compressed --memory-mapping \
		--report $@ --output /dev/null \
		$(TRIMMED_DIR)/$*_1.fastq.gz $(TRIMMED_DIR)/$*_2.fastq.gz

# QUAST ------------------------------------------------------------------------
$(QUAST_DIR)/.done: $(patsubst %,$(QUAST_DIR)/%/report.tsv,$(SAMPLES)) | $(QUAST_DIR)
	@touch $@

$(QUAST_DIR)/%/report.tsv: $(ASSEMBLY_DIR)/%_assembled.fasta | $(QUAST_DIR)
	$(RUN_BIOINFO) quast.py \
		--output-dir $(QUAST_DIR)/$* \
		--threads $(QUAST_THREADS) \
		$<

$(TRIMMED_DIR) $(QC_DIR) $(KRAKEN_DIR) $(QUAST_DIR) $(CLEAN_DIR):
	mkdir -p $@

# Assembly ----------------------------------------------------------------------
$(DOWNSAMPLED_DIR)/.done: $(patsubst %,$(DOWNSAMPLED_DIR)/%_1.fastq.gz,$(SAMPLES)) | $(DOWNSAMPLED_DIR)
	@touch $@

$(DOWNSAMPLED_DIR)/%_1.fastq.gz $(DOWNSAMPLED_DIR)/%_2.fastq.gz &: $(TRIMMED_DIR)/%_1.fastq.gz $(TRIMMED_DIR)/%_2.fastq.gz $(QC_DIR)/%_fastp.json | $(DOWNSAMPLED_DIR)
	fraction=$$($(RUN_AMR) python3 -c "import json; d=json.load(open('$(QC_DIR)/$*_fastp.json')); cov=d['summary']['after_filtering']['total_bases']/$(GENOME_SIZE); print(min(0.999999, $(TARGET_COVERAGE)/cov))"); \
	$(RUN_BIOINFO) seqtk sample -s42 $(TRIMMED_DIR)/$*_1.fastq.gz $$fraction | gzip > $(DOWNSAMPLED_DIR)/$*_1.fastq.gz; \
	$(RUN_BIOINFO) seqtk sample -s42 $(TRIMMED_DIR)/$*_2.fastq.gz $$fraction | gzip > $(DOWNSAMPLED_DIR)/$*_2.fastq.gz

$(ASSEMBLY_DIR)/.done: $(patsubst %,$(ASSEMBLY_DIR)/%_assembled.fasta,$(SAMPLES)) | $(ASSEMBLY_DIR)
	@touch $@

$(ASSEMBLY_DIR)/%_assembled.fasta: $(DOWNSAMPLED_DIR)/%_1.fastq.gz $(DOWNSAMPLED_DIR)/%_2.fastq.gz | $(ASSEMBLY_DIR)
	$(RUN_BIOINFO) spades.py \
		--pe1-1 $(DOWNSAMPLED_DIR)/$*_1.fastq.gz \
		--pe1-2 $(DOWNSAMPLED_DIR)/$*_2.fastq.gz \
		--isolate --only-assembler \
		-k 21,33,55 \
		--threads $(SPADES_THREADS) --memory $(SPADES_MEMORY) \
		-o $(ASSEMBLY_DIR)/$*_tmp
	mv $(ASSEMBLY_DIR)/$*_tmp/contigs.fasta $@
	rm -rf $(ASSEMBLY_DIR)/$*_tmp

$(DOWNSAMPLED_DIR) $(ASSEMBLY_DIR):
	mkdir -p $@

# AMRFinderPlus ----------------------------------------------------------------
$(AMR_DIR)/.done: $(patsubst %,$(AMR_DIR)/%_amr.tsv,$(SAMPLES)) | $(AMR_DIR)
	@touch $@

$(AMRFINDER_DB)/latest/AMR.LIB:
	mkdir -p $(AMRFINDER_DB)
	$(RUN_BIOINFO) amrfinder_update -d $(AMRFINDER_DB)

$(AMR_DIR)/%_amr.tsv $(AMR_DIR)/%_amr_genes.fna &: $(ASSEMBLY_DIR)/%_assembled.fasta $(AMRFINDER_DB)/latest/AMR.LIB | $(AMR_DIR)
	$(RUN_BIOINFO) amrfinder \
		--nucleotide $< \
		--organism Klebsiella_pneumoniae \
		--database $(AMRFINDER_DB)/latest \
		--output $(AMR_DIR)/$*_amr.tsv \
		--nucleotide_output $(AMR_DIR)/$*_amr_genes.fna \
		--threads $(AMRFINDER_THREADS)
	touch $(AMR_DIR)/$*_amr_genes.fna

$(AMR_DIR):
	mkdir -p $@

# Per-sample cleanup of large intermediates ------------------------------------
# Once AMR, Kraken2, and QUAST are done for a sample, delete its big files
# (raw reads, trimmed reads, downsampled reads, assembly) before moving on.
# AMR/QC outputs are intentionally kept here; they are removed later once
# features and the aggregated MultiQC report are built.
$(CLEAN_DIR)/%_cleaned: $(AMR_DIR)/%_amr.tsv $(KRAKEN_DIR)/%_report.txt $(QUAST_DIR)/%/report.tsv | $(CLEAN_DIR)
	rm -f $(READS_DIR)/$*_1.fastq.gz $(READS_DIR)/$*_2.fastq.gz \
	      $(TRIMMED_DIR)/$*_1.fastq.gz $(TRIMMED_DIR)/$*_2.fastq.gz \
	      $(DOWNSAMPLED_DIR)/$*_1.fastq.gz $(DOWNSAMPLED_DIR)/$*_2.fastq.gz \
	      $(ASSEMBLY_DIR)/$*_assembled.fasta
	@touch $@

$(QC_DIR)/.done: $(patsubst %,$(QC_DIR)/%_fastp.json,$(SAMPLES)) | $(QC_DIR)
	@touch $@

# Features (tabular) -----------------------------------------------------------
$(FEATURES_DIR)/train_features.csv $(FEATURES_DIR)/test_features.csv: $(FEATURES_DIR)/.done
	@true

$(FEATURES_DIR)/.done: $(patsubst %,$(CLEAN_DIR)/%_cleaned,$(SAMPLES)) results/metadata/train.csv results/metadata/test.csv scripts/build_features.py | $(FEATURES_DIR)
	$(RUN_AMR) python3 scripts/build_features.py \
		--amr-files $(AMR_DIR)/*_amr.tsv \
		--train-metadata results/metadata/train.csv \
		--test-metadata results/metadata/test.csv \
		--train-output $(FEATURES_DIR)/train_features.csv \
		--test-output $(FEATURES_DIR)/test_features.csv
	touch $@

$(FEATURES_DIR):
	mkdir -p $@

# Sequences (for optional DNABERT-2) -------------------------------------------
$(SEQUENCES_DIR)/train_sequences.csv $(SEQUENCES_DIR)/test_sequences.csv: $(SEQUENCES_DIR)/.done
	@true

$(SEQUENCES_DIR)/.done: $(patsubst %,$(CLEAN_DIR)/%_cleaned,$(SAMPLES)) results/metadata/train.csv results/metadata/test.csv scripts/build_sequences.py | $(SEQUENCES_DIR)
	$(RUN_AMR) python3 scripts/build_sequences.py \
		--seq-files $(AMR_DIR)/*_amr_genes.fna \
		--train-metadata results/metadata/train.csv \
		--test-metadata results/metadata/test.csv \
		--train-output $(SEQUENCES_DIR)/train_sequences.csv \
		--test-output $(SEQUENCES_DIR)/test_sequences.csv
	touch $@

$(SEQUENCES_DIR):
	mkdir -p $@

# Hyperparameter tuning (Optuna) -----------------------------------------------
tune: $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/xgboost/tune_$(abx).json) \
      $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/lightgbm/tune_$(abx).json) \
      $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/nn/tune_$(abx).json)

$(MODELS_DIR)/xgboost/tune_%.json: $(FEATURES_DIR)/train_features.csv scripts/tune_xgboost.py | $(MODELS_DIR)/xgboost
	$(RUN_AMR) python3 scripts/tune_xgboost.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--n-trials $(TUNE_TRIALS_XGB) \
		--n-splits $(TUNE_SPLITS) \
		--output $@

$(MODELS_DIR)/lightgbm/tune_%.json: $(FEATURES_DIR)/train_features.csv scripts/tune_lightgbm.py | $(MODELS_DIR)/lightgbm
	$(RUN_AMR) python3 scripts/tune_lightgbm.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--n-trials $(TUNE_TRIALS_LGB) \
		--n-splits $(TUNE_SPLITS) \
		--output $@

$(MODELS_DIR)/nn/tune_%.json: $(FEATURES_DIR)/train_features.csv scripts/tune_nn.py | $(MODELS_DIR)/nn
	$(RUN_AMR) python3 scripts/tune_nn.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--n-trials $(TUNE_TRIALS_NN) \
		--n-splits $(TUNE_SPLITS) \
		--output $@

# Models -----------------------------------------------------------------------
models:
	@if [ -z "$(SAMPLES)" ]; then \
		$(MAKE) metadata MAX_SAMPLES=$(MAX_SAMPLES) && $(MAKE) models MAX_SAMPLES=$(MAX_SAMPLES); \
	else \
		$(MAKE) $(MODELS_DIR)/.done; \
	fi

$(MODELS_DIR)/.done: $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/xgboost/$(abx)_metrics.json) \
                     $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/lightgbm/$(abx)_metrics.json) \
                     $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/nn/$(abx)_metrics.json)
	@touch $@

$(MODELS_DIR)/xgboost/%_metrics.json: $(FEATURES_DIR)/train_features.csv $(FEATURES_DIR)/test_features.csv scripts/train_xgboost.py $(MODELS_DIR)/xgboost/tune_%.json | $(MODELS_DIR)/xgboost
	$(RUN_AMR) python3 scripts/train_xgboost.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--test-features $(FEATURES_DIR)/test_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--params-input $(MODELS_DIR)/xgboost/tune_$*.json \
		--model-output $(MODELS_DIR)/xgboost/$*_model.json \
		--params-output $(MODELS_DIR)/xgboost/$*_params.json \
		--metrics-output $@ \
		--predictions-output $(MODELS_DIR)/xgboost/$*_predictions.csv \
		--importance-output $(MODELS_DIR)/xgboost/$*_importance.csv \
		--shap-plot-output $(MODELS_DIR)/xgboost/$*_shap.png

$(MODELS_DIR)/lightgbm/%_metrics.json: $(FEATURES_DIR)/train_features.csv $(FEATURES_DIR)/test_features.csv scripts/train_lightgbm.py $(MODELS_DIR)/lightgbm/tune_%.json | $(MODELS_DIR)/lightgbm
	$(RUN_AMR) python3 scripts/train_lightgbm.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--test-features $(FEATURES_DIR)/test_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--params-input $(MODELS_DIR)/lightgbm/tune_$*.json \
		--model-output $(MODELS_DIR)/lightgbm/$*_model.txt \
		--params-output $(MODELS_DIR)/lightgbm/$*_params.json \
		--metrics-output $@ \
		--predictions-output $(MODELS_DIR)/lightgbm/$*_predictions.csv \
		--importance-output $(MODELS_DIR)/lightgbm/$*_importance.csv \
		--shap-plot-output $(MODELS_DIR)/lightgbm/$*_shap.png

$(MODELS_DIR)/nn/%_metrics.json: $(FEATURES_DIR)/train_features.csv $(FEATURES_DIR)/test_features.csv scripts/train_nn.py $(MODELS_DIR)/nn/tune_%.json | $(MODELS_DIR)/nn
	$(RUN_AMR) python3 scripts/train_nn.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--test-features $(FEATURES_DIR)/test_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--params-input $(MODELS_DIR)/nn/tune_$*.json \
		--model-output $(MODELS_DIR)/nn/$*_model.pt \
		--params-output $(MODELS_DIR)/nn/$*_params.json \
		--metrics-output $@ \
		--predictions-output $(MODELS_DIR)/nn/$*_predictions.csv \
		--importance-output $(MODELS_DIR)/nn/$*_importance.csv \
		--importance-plot-output $(MODELS_DIR)/nn/$*_importance.png

$(MODELS_DIR)/xgboost $(MODELS_DIR)/lightgbm $(MODELS_DIR)/nn:
	mkdir -p $@

# DNABERT-2 (optional) ---------------------------------------------------------
dnabert:
	@if [ -z "$(SAMPLES)" ]; then \
		$(MAKE) metadata MAX_SAMPLES=$(MAX_SAMPLES) && $(MAKE) dnabert MAX_SAMPLES=$(MAX_SAMPLES); \
	else \
		$(MAKE) $(MODELS_DIR)/dnabert/.done; \
	fi

$(MODELS_DIR)/dnabert/.done: $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/dnabert/$(abx)_metrics.json) | $(MODELS_DIR)/dnabert
	@touch $@

$(MODELS_DIR)/dnabert/%_metrics.json: $(SEQUENCES_DIR)/train_sequences.csv $(SEQUENCES_DIR)/test_sequences.csv scripts/train_dnabert.py | $(MODELS_DIR)/dnabert
	$(RUN_AMR) python3 scripts/train_dnabert.py \
		--train-features $(SEQUENCES_DIR)/train_sequences.csv \
		--test-features $(SEQUENCES_DIR)/test_sequences.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--model-output $(MODELS_DIR)/dnabert/$*_model.pt \
		--params-output $(MODELS_DIR)/dnabert/$*_params.json \
		--metrics-output $@ \
		--predictions-output $(MODELS_DIR)/dnabert/$*_predictions.csv \
		--importance-output $(MODELS_DIR)/dnabert/$*_importance.csv \
		--importance-plot-output $(MODELS_DIR)/dnabert/$*_importance.png

$(MODELS_DIR)/dnabert:
	mkdir -p $@

# MultiQC: aggregate fastp + kraken2 + quast -----------------------------------
multiqc:
	@if [ -z "$(SAMPLES)" ]; then \
		$(MAKE) metadata MAX_SAMPLES=$(MAX_SAMPLES) && $(MAKE) multiqc MAX_SAMPLES=$(MAX_SAMPLES); \
	else \
		$(MAKE) $(MULTIQC_DIR)/multiqc_report.html; \
	fi

$(MULTIQC_DIR)/multiqc_report.html: $(QC_DIR)/.done $(KRAKEN_DIR)/.done $(QUAST_DIR)/.done | $(MULTIQC_DIR)
	$(RUN_BIOINFO) multiqc $(QC_DIR) $(KRAKEN_DIR) $(QUAST_DIR) \
		--outdir $(MULTIQC_DIR) \
		--filename multiqc_report.html \
		--force

$(MULTIQC_DIR):
	mkdir -p $@

# Report / interpretability summary --------------------------------------------
report:
	@if [ -z "$(SAMPLES)" ]; then \
		$(MAKE) metadata MAX_SAMPLES=$(MAX_SAMPLES) && $(MAKE) report MAX_SAMPLES=$(MAX_SAMPLES); \
	else \
		$(MAKE) $(REPORT_DIR)/summary.json $(MULTIQC_DIR)/multiqc_report.html && \
		rm -rf $(DATA_DIR) $(FEATURES_DIR) $(SEQUENCES_DIR); \
	fi

$(REPORT_DIR)/summary.json: $(MODELS_DIR)/.done scripts/summarize.py | $(REPORT_DIR)
	$(RUN_AMR) python3 scripts/summarize.py \
		--models-dir $(MODELS_DIR) \
		--antibiotics $(ANTIBIOTICS) \
		--output $@

$(REPORT_DIR):
	mkdir -p $@

# Clean ------------------------------------------------------------------------
clean:
	rm -rf $(DATA_DIR) $(RESULTS_DIR) reference
