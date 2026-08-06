# Klebsiella pneumoniae AMR prediction pipeline
# Orchestrated with GNU Make + conda.

# macOS ships BSD make by default; this Makefile uses GNU Make features
# (group targets, pattern rules). On macOS install gmake and run 'gmake'.
ifeq ($(findstring .,$(MAKE_VERSION)),)
  $(error This Makefile requires GNU Make. On macOS: brew install make && gmake)
endif

# Preserve intermediate files (reads, trimmed, assemblies, etc.)
.SECONDARY:

# Auto-generated includes ------------------------------------------------------
# Only include generated config if the amr environment exists; otherwise
# 'make setup' would try to build config.mk before the environment exists.
CONDA_BASE := $(shell conda info --base 2>/dev/null)
AMR_ENV := $(wildcard $(CONDA_BASE)/envs/amr)
ifneq ($(AMR_ENV),)
  -include results/metadata/config.mk
  -include results/metadata/samples.mk
endif

# Conda runners ----------------------------------------------------------------
RUN_AMR     := conda run --no-capture-output -n amr
RUN_BIOINFO := conda run --no-capture-output -n bioinfo

# Development sample cap. Override per run, e.g.:
#   make metadata MAX_SAMPLES=5
#   make metadata MAX_SAMPLES=20
#   make metadata MAX_SAMPLES=-1   # full dataset
MAX_SAMPLES := -1

# Directories ------------------------------------------------------------------
RESULTS_DIR        := results
DATA_DIR           := data
READS_DIR          := $(DATA_DIR)/reads
TRIMMED_DIR        := $(DATA_DIR)/trimmed
QC_DIR             := $(DATA_DIR)/qc
KRAKEN_DIR         := $(DATA_DIR)/kraken2
DOWNSAMPLED_DIR    := $(DATA_DIR)/downsampled
ASSEMBLY_DIR       := $(DATA_DIR)/assembly
ASSEMBLY_STATS_DIR := $(DATA_DIR)/assembly_stats
AMR_DIR            := $(DATA_DIR)/amr
MULTIQC_DIR        := $(RESULTS_DIR)/multiqc
FEATURES_DIR       := $(RESULTS_DIR)/features
SEQUENCES_DIR      := $(RESULTS_DIR)/sequences
MODELS_DIR         := $(RESULTS_DIR)/models
REPORT_DIR         := $(RESULTS_DIR)/report

# Phony targets ----------------------------------------------------------------
.PHONY: all setup test metadata reads qc kraken2 assembly amr features sequences models dnabert report clean

all:
	@if [ -z "$(SAMPLES)" ]; then \
		echo "ERROR: run 'make metadata MAX_SAMPLES=...' first"; exit 1; \
	fi
	$(MAKE) report

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
reads: $(patsubst %,$(READS_DIR)/%_1.fastq.gz,$(SAMPLES))

$(READS_DIR)/%_1.fastq.gz $(READS_DIR)/%_2.fastq.gz &: scripts/download_reads.py | $(READS_DIR)
	$(RUN_AMR) python3 scripts/download_reads.py \
		--accession $* \
		--out1 $(READS_DIR)/$*_1.fastq.gz \
		--out2 $(READS_DIR)/$*_2.fastq.gz

$(READS_DIR):
	mkdir -p $@

# QC (fastp + MultiQC) ---------------------------------------------------------
qc: $(patsubst %,$(QC_DIR)/%_fastp.json,$(SAMPLES)) $(MULTIQC_DIR)/fastp_multiqc_report.html

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

$(MULTIQC_DIR)/fastp_multiqc_report.html: $(patsubst %,$(QC_DIR)/%_fastp.json,$(SAMPLES)) | $(MULTIQC_DIR)
	$(RUN_BIOINFO) multiqc $(QC_DIR) --outdir $(MULTIQC_DIR) --filename fastp_multiqc_report.html --force

$(TRIMMED_DIR) $(QC_DIR) $(MULTIQC_DIR):
	mkdir -p $@

# Kraken2 contamination screening (optional, requires ~8 GB database) ----------
kraken2: $(patsubst %,$(KRAKEN_DIR)/%_report.txt,$(SAMPLES))

$(KRAKEN_DB)/taxo.k2d:
	mkdir -p $(KRAKEN_DB)
	curl -L -o $(KRAKEN_DB)/k2_standard_08gb_20250402.tar.gz \
		https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz
	tar -xzf $(KRAKEN_DB)/k2_standard_08gb_20250402.tar.gz -C $(KRAKEN_DB)
	rm $(KRAKEN_DB)/k2_standard_08gb_20250402.tar.gz

$(KRAKEN_DIR)/%_report.txt: $(TRIMMED_DIR)/%_1.fastq.gz $(TRIMMED_DIR)/%_2.fastq.gz $(KRAKEN_DB)/taxo.k2d | $(KRAKEN_DIR)
	$(RUN_BIOINFO) kraken2 --db $(KRAKEN_DB) \
		--paired --gzip-compressed \
		--report $@ --output /dev/null \
		$(TRIMMED_DIR)/$*_1.fastq.gz $(TRIMMED_DIR)/$*_2.fastq.gz

$(KRAKEN_DIR):
	mkdir -p $@

# Assembly ----------------------------------------------------------------------
assembly: $(patsubst %,$(ASSEMBLY_DIR)/%_assembled.fasta,$(SAMPLES)) \
          $(patsubst %,$(ASSEMBLY_STATS_DIR)/%_stats.tsv,$(SAMPLES))

$(DOWNSAMPLED_DIR)/%_1.fastq.gz $(DOWNSAMPLED_DIR)/%_2.fastq.gz &: $(TRIMMED_DIR)/%_1.fastq.gz $(TRIMMED_DIR)/%_2.fastq.gz $(QC_DIR)/%_fastp.json | $(DOWNSAMPLED_DIR)
	fraction=$$($(RUN_AMR) python3 -c "import json; d=json.load(open('$(QC_DIR)/$*_fastp.json')); cov=d['summary']['after_filtering']['total_bases']/$(GENOME_SIZE); print(min(0.999999, $(TARGET_COVERAGE)/cov))"); \
	$(RUN_BIOINFO) seqtk sample -s42 $(TRIMMED_DIR)/$*_1.fastq.gz $$fraction | gzip > $(DOWNSAMPLED_DIR)/$*_1.fastq.gz; \
	$(RUN_BIOINFO) seqtk sample -s42 $(TRIMMED_DIR)/$*_2.fastq.gz $$fraction | gzip > $(DOWNSAMPLED_DIR)/$*_2.fastq.gz

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

$(ASSEMBLY_STATS_DIR)/%_stats.tsv: $(ASSEMBLY_DIR)/%_assembled.fasta | $(ASSEMBLY_STATS_DIR)
	$(RUN_AMR) python3 scripts/assembly_stats.py --assembly $< --output $@

$(DOWNSAMPLED_DIR) $(ASSEMBLY_DIR) $(ASSEMBLY_STATS_DIR):
	mkdir -p $@

# AMRFinderPlus ----------------------------------------------------------------
amr: $(patsubst %,$(AMR_DIR)/%_amr.tsv,$(SAMPLES))

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

# Features (tabular) -----------------------------------------------------------
features: $(FEATURES_DIR)/.done

$(FEATURES_DIR)/train_features.csv $(FEATURES_DIR)/test_features.csv: $(FEATURES_DIR)/.done
	@true

$(FEATURES_DIR)/.done: $(patsubst %,$(AMR_DIR)/%_amr.tsv,$(SAMPLES)) results/metadata/train.csv results/metadata/test.csv scripts/build_features.py | $(FEATURES_DIR)
	$(RUN_AMR) python3 scripts/build_features.py \
		--amr-files $(wildcard $(AMR_DIR)/*_amr.tsv) \
		--train-metadata results/metadata/train.csv \
		--test-metadata results/metadata/test.csv \
		--train-output $(FEATURES_DIR)/train_features.csv \
		--test-output $(FEATURES_DIR)/test_features.csv
	touch $@

$(FEATURES_DIR):
	mkdir -p $@

# Sequences (for optional DNABERT-2) -------------------------------------------
sequences: $(SEQUENCES_DIR)/.done

$(SEQUENCES_DIR)/train_sequences.csv $(SEQUENCES_DIR)/test_sequences.csv: $(SEQUENCES_DIR)/.done
	@true

$(SEQUENCES_DIR)/.done: $(patsubst %,$(AMR_DIR)/%_amr_genes.fna,$(SAMPLES)) results/metadata/train.csv results/metadata/test.csv scripts/build_sequences.py | $(SEQUENCES_DIR)
	$(RUN_AMR) python3 scripts/build_sequences.py \
		--seq-files $(wildcard $(AMR_DIR)/*_amr_genes.fna) \
		--train-metadata results/metadata/train.csv \
		--test-metadata results/metadata/test.csv \
		--train-output $(SEQUENCES_DIR)/train_sequences.csv \
		--test-output $(SEQUENCES_DIR)/test_sequences.csv
	touch $@

$(SEQUENCES_DIR):
	mkdir -p $@

# Models -----------------------------------------------------------------------
models: $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/xgboost/$(abx)_metrics.json) \
        $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/lightgbm/$(abx)_metrics.json) \
        $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/nn/$(abx)_metrics.json)

$(MODELS_DIR)/xgboost/%_metrics.json: $(FEATURES_DIR)/train_features.csv $(FEATURES_DIR)/test_features.csv scripts/train_xgboost.py | $(MODELS_DIR)/xgboost
	$(RUN_AMR) python3 scripts/train_xgboost.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--test-features $(FEATURES_DIR)/test_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--model-output $(MODELS_DIR)/xgboost/$*_model.json \
		--params-output $(MODELS_DIR)/xgboost/$*_params.json \
		--metrics-output $@ \
		--predictions-output $(MODELS_DIR)/xgboost/$*_predictions.csv \
		--importance-output $(MODELS_DIR)/xgboost/$*_importance.csv \
		--shap-plot-output $(MODELS_DIR)/xgboost/$*_shap.png

$(MODELS_DIR)/lightgbm/%_metrics.json: $(FEATURES_DIR)/train_features.csv $(FEATURES_DIR)/test_features.csv scripts/train_lightgbm.py | $(MODELS_DIR)/lightgbm
	$(RUN_AMR) python3 scripts/train_lightgbm.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--test-features $(FEATURES_DIR)/test_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--model-output $(MODELS_DIR)/lightgbm/$*_model.txt \
		--params-output $(MODELS_DIR)/lightgbm/$*_params.json \
		--metrics-output $@ \
		--predictions-output $(MODELS_DIR)/lightgbm/$*_predictions.csv \
		--importance-output $(MODELS_DIR)/lightgbm/$*_importance.csv \
		--shap-plot-output $(MODELS_DIR)/lightgbm/$*_shap.png

$(MODELS_DIR)/nn/%_metrics.json: $(FEATURES_DIR)/train_features.csv $(FEATURES_DIR)/test_features.csv scripts/train_nn.py | $(MODELS_DIR)/nn
	$(RUN_AMR) python3 scripts/train_nn.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--test-features $(FEATURES_DIR)/test_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--model-output $(MODELS_DIR)/nn/$*_model.pt \
		--params-output $(MODELS_DIR)/nn/$*_params.json \
		--metrics-output $@ \
		--predictions-output $(MODELS_DIR)/nn/$*_predictions.csv \
		--importance-output $(MODELS_DIR)/nn/$*_importance.csv \
		--importance-plot-output $(MODELS_DIR)/nn/$*_importance.png

$(MODELS_DIR)/xgboost $(MODELS_DIR)/lightgbm $(MODELS_DIR)/nn:
	mkdir -p $@

# DNABERT-2 (optional) ---------------------------------------------------------
dnabert: $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/dnabert/$(abx)_metrics.json)

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

# Report / interpretability summary --------------------------------------------
report: $(REPORT_DIR)/summary.json

$(REPORT_DIR)/summary.json: models scripts/summarize.py | $(REPORT_DIR)
	$(RUN_AMR) python3 scripts/summarize.py \
		--models-dir $(MODELS_DIR) \
		--antibiotics $(ANTIBIOTICS) \
		--output $@

$(REPORT_DIR):
	mkdir -p $@

# Clean ------------------------------------------------------------------------
clean:
	rm -rf $(DATA_DIR) $(RESULTS_DIR) reference
