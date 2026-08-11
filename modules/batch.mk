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
