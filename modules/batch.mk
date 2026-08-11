# Batch processing (disk-efficient) -------------------------------------------
# Downloads all samples in parallel (skipping inaccessible ones), then processes
# each sample end-to-end through all stages (trim -> downsample -> assemble ->
# amr -> kraken2 -> quast) before cleaning up intermediate files. This limits
# peak disk usage to ~BATCH_SIZE samples' worth of intermediate data.

_download-all: $(READS_DIR)
	@echo "Downloading reads for $(words $(SAMPLES)) samples ($(DOWNLOAD_JOBS) at a time)..."
	@echo "$(SAMPLES)" | tr ' ' '\n' | \
	xargs -P $(DOWNLOAD_JOBS) -I {} sh -c '$(MAKE) --no-print-directory $(READS_DIR)/{}_1.fastq.gz 2>&1 || true'
	@skipped=$$(ls $(READS_DIR)/*.skip 2>/dev/null | wc -l); \
	active=$$(ls $(READS_DIR)/*_1.fastq.gz 2>/dev/null | wc -l); \
	echo "Downloaded: $$active samples | Skipped: $$skipped samples"

_process-samples:
	@if [ -z "$(SAMPLES_WITH_READS)" ]; then \
		echo "No samples with reads to process."; \
		exit 0; \
	fi
	@echo "Processing $(words $(SAMPLES_WITH_READS)) samples end-to-end ($(BATCH_SIZE) at a time)..."
	@echo "$(SAMPLES_WITH_READS)" | tr ' ' '\n' | \
	xargs -P $(BATCH_SIZE) -I {} sh -c ' \
		echo "=== {} ==="; \
		$(MAKE) --no-print-directory $(CLEAN_DIR)/{}_cleaned && echo "{}: OK" || \
		(echo "{}: FAILED"; rm -f $(READS_DIR)/{}_1.fastq.gz $(READS_DIR)/{}_2.fastq.gz; mkdir -p $(CLEAN_DIR); touch $(CLEAN_DIR)/{}_cleaned); \
		true; \
	'

process-samples:
	@if [ -z "$(SAMPLES)" ]; then \
		$(MAKE) --no-print-directory metadata MAX_SAMPLES=$(MAX_SAMPLES) && \
		$(MAKE) --no-print-directory process-samples MAX_SAMPLES=$(MAX_SAMPLES); \
	else \
		$(MAKE) --no-print-directory _download-all && \
		$(MAKE) --no-print-directory _process-samples; \
	fi

# MultiQC: aggregate fastp + kraken2 + quast ----------------------------------
multiqc:
	@if [ -z "$(SAMPLES)" ]; then \
		$(MAKE) --no-print-directory metadata MAX_SAMPLES=$(MAX_SAMPLES) && $(MAKE) --no-print-directory multiqc MAX_SAMPLES=$(MAX_SAMPLES); \
	else \
		$(MAKE) --no-print-directory $(MULTIQC_DIR)/multiqc_report.html; \
	fi

$(MULTIQC_DIR)/multiqc_report.html: $(QC_DIR)/.done $(KRAKEN_DIR)/.done $(QUAST_DIR)/.done | $(MULTIQC_DIR)
	$(RUN_BIOINFO) multiqc $(QC_DIR) $(KRAKEN_DIR) $(QUAST_DIR) \
		--outdir $(MULTIQC_DIR) \
		--filename multiqc_report.html \
		--force

# Post-assembly analysis (features + sequences + models + multiqc + summary) ---
analyze:
	@if [ -z "$(SAMPLES)" ]; then \
		$(MAKE) --no-print-directory metadata MAX_SAMPLES=$(MAX_SAMPLES) && \
		$(MAKE) --no-print-directory analyze MAX_SAMPLES=$(MAX_SAMPLES); \
	else \
		$(MAKE) --no-print-directory $(FEATURES_DIR)/.done \
		                   $(SEQUENCES_DIR)/.done \
		                   $(MODELS_DIR)/.done \
		                   $(MULTIQC_DIR)/multiqc_report.html \
		                   $(REPORT_DIR)/summary.json; \
	fi

# Report / interpretability summary -------------------------------------------
$(REPORT_DIR)/summary.json: $(MODELS_DIR)/.done scripts/summarize.py | $(REPORT_DIR)
	$(RUN_AMR) python3 scripts/summarize.py \
		--models-dir $(MODELS_DIR) \
		--antibiotics $(ANTIBIOTICS) \
		--output $@

report:
	@if [ -z "$(SAMPLES)" ]; then \
		$(MAKE) --no-print-directory metadata MAX_SAMPLES=$(MAX_SAMPLES) && $(MAKE) --no-print-directory report MAX_SAMPLES=$(MAX_SAMPLES); \
	else \
		$(MAKE) --no-print-directory process-samples && \
		$(MAKE) --no-print-directory analyze && \
		rm -rf $(READS_DIR) $(TRIMMED_DIR) $(DOWNSAMPLED_DIR) $(QC_DIR) $(KRAKEN_DIR) $(QUAST_DIR) $(AMR_DIR) $(CLEAN_DIR) $(FEATURES_DIR) $(SEQUENCES_DIR); \
	fi
