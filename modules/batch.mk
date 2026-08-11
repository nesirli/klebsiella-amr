# Batch processing (disk-efficient) -------------------------------------------
# Each sample is downloaded and processed end-to-end (trim -> downsample ->
# assemble -> amr -> kraken2 -> quast -> cleanup) within a single batch.
# BATCH_SIZE controls both download and processing concurrency. At most
# BATCH_SIZE samples' worth of intermediate data ever exists on disk.

_process-samples:
	@echo "Processing $(words $(SAMPLES)) samples ($(BATCH_SIZE) at a time)..."
	@batch=0; \
	for sample in $(SAMPLES); do \
		( \
			sample=$$sample; \
			if [ -f $(READS_DIR)/$${sample}.skip ]; then \
				mkdir -p $(CLEAN_DIR); \
				touch $(CLEAN_DIR)/$${sample}_cleaned; \
				echo "$$sample: SKIPPED"; \
			else \
				{ [ -f $(READS_DIR)/$${sample}_1.fastq.gz ] || \
				  $(MAKE) --no-print-directory $(READS_DIR)/$${sample}_1.fastq.gz 2>&1; } || true; \
				if [ -f $(READS_DIR)/$${sample}.skip ]; then \
					mkdir -p $(CLEAN_DIR); \
					touch $(CLEAN_DIR)/$${sample}_cleaned; \
					echo "$$sample: SKIPPED"; \
				elif [ -f $(READS_DIR)/$${sample}_1.fastq.gz ]; then \
					echo "=== $$sample ==="; \
					if $(MAKE) --no-print-directory $(CLEAN_DIR)/$${sample}_cleaned; then \
						echo "$$sample: OK"; \
					else \
						echo "$$sample: FAILED"; \
						rm -f $(READS_DIR)/$${sample}_1.fastq.gz $(READS_DIR)/$${sample}_2.fastq.gz; \
						mkdir -p $(CLEAN_DIR); \
						touch $(CLEAN_DIR)/$${sample}_cleaned; \
					fi; \
				else \
					mkdir -p $(CLEAN_DIR); \
					touch $(CLEAN_DIR)/$${sample}_cleaned; \
					echo "$$sample: DOWNLOAD FAILED"; \
				fi; \
			fi \
		) & \
		batch=$$((batch + 1)); \
		if [ $$batch -ge $(BATCH_SIZE) ]; then wait; batch=0; fi; \
	done; \
	wait; \
	true

process-samples:
	@if [ -z "$(SAMPLES)" ]; then \
		$(MAKE) --no-print-directory metadata MAX_SAMPLES=$(MAX_SAMPLES) && \
		$(MAKE) --no-print-directory process-samples MAX_SAMPLES=$(MAX_SAMPLES); \
	else \
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
