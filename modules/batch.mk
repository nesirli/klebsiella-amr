# Batch processing (disk-efficient) -------------------------------------------
# Each sample is downloaded and processed end-to-end (trim -> downsample ->
# assemble -> amr -> kraken2 -> quast -> cleanup) within a single batch.
# BATCH_SIZE controls both download and processing concurrency. At most
# BATCH_SIZE samples' worth of intermediate data ever exists on disk.

# Retiring a sample removes every intermediate it owns, INCLUDING its fastp
# JSON. SAMPLES_ACTIVE is derived from those JSONs, so a sample that could not
# be downloaded or processed drops out of the run entirely. Leaving the JSON
# behind would keep the sample in SAMPLES_ACTIVE, and the analyze phase would
# then demand its quast report and walk all the way back to re-downloading
# reads this recipe had just deleted -- failing the whole run over one bad
# sample. Assemblies are deliberately kept elsewhere, but a half-written one
# from a crashed SPAdes is not worth keeping.
define retire-sample
rm -rf $(QUAST_DIR)/$$s $(ASSEMBLY_DIR)/$${s}_tmp; \
rm -f $(READS_DIR)/$${s}_1.fastq.gz $(READS_DIR)/$${s}_2.fastq.gz \
      $(TRIMMED_DIR)/$${s}_1.fastq.gz $(TRIMMED_DIR)/$${s}_2.fastq.gz \
      $(DOWNSAMPLED_DIR)/$${s}_1.fastq.gz $(DOWNSAMPLED_DIR)/$${s}_2.fastq.gz \
      $(QC_DIR)/$${s}_fastp.json $(QC_DIR)/$${s}_fastp.html \
      $(KRAKEN_DIR)/$${s}_report.txt \
      $(AMR_DIR)/$${s}_amr.tsv $(AMR_DIR)/$${s}_amr_genes.fna
endef

_process-samples:
	@echo "Processing $(words $(SAMPLES)) samples ($(BATCH_SIZE) at a time)..."
	@batch=0; \
	for sample in $(SAMPLES); do \
		$(MAKE) --no-print-directory _process-one SAMPLE=$$sample & \
		batch=$$((batch + 1)); \
		if [ $$batch -ge $(BATCH_SIZE) ]; then wait; batch=0; fi; \
	done; \
	wait

# One sample, download through cleanup. Always exits 0: a sample that cannot be
# fetched or processed is retired rather than failing the batch, so one bad
# accession out of thousands does not sink the run.
_process-one:
	@s='$(SAMPLE)'; \
	if [ -f $(READS_DIR)/$$s.skip ]; then \
		echo "$$s: SKIPPED"; exit 0; \
	fi; \
	if [ ! -f $(READS_DIR)/$${s}_1.fastq.gz ]; then \
		$(MAKE) --no-print-directory $(READS_DIR)/$${s}_1.fastq.gz 2>&1 || true; \
	fi; \
	if [ -f $(READS_DIR)/$$s.skip ]; then \
		echo "$$s: SKIPPED"; exit 0; \
	fi; \
	if [ ! -f $(READS_DIR)/$${s}_1.fastq.gz ]; then \
		echo "$$s: DOWNLOAD FAILED"; $(retire-sample); exit 0; \
	fi; \
	echo "=== $$s ==="; \
	if $(MAKE) --no-print-directory $(CLEAN_DIR)/$${s}_cleaned; then \
		echo "$$s: OK"; \
	else \
		echo "$$s: FAILED"; $(retire-sample); \
	fi

process-samples:
	$(call need-samples,_process-samples)

# MultiQC: aggregate fastp + kraken2 + quast ----------------------------------
multiqc:
	$(call need-samples,$(MULTIQC_DIR)/multiqc_report.html)

$(MULTIQC_DIR)/multiqc_report.html: $(QC_DIR)/.done $(KRAKEN_DIR)/.done $(QUAST_DIR)/.done | $(MULTIQC_DIR)
	$(RUN_BIOINFO) multiqc $(QC_DIR) $(KRAKEN_DIR) $(QUAST_DIR) \
		--outdir $(MULTIQC_DIR) \
		--filename multiqc_report.html \
		--force

# Post-assembly analysis (features + sequences + models + multiqc + summary) ---
analyze:
	$(call need-samples,$(FEATURES_DIR)/train_features.csv \
	                    $(SEQUENCES_DIR)/train_sequences.csv \
	                    $(MODELS_DIR)/.done \
	                    $(MULTIQC_DIR)/multiqc_report.html \
	                    $(REPORT_DIR)/summary.json)

# Report / interpretability summary -------------------------------------------
$(REPORT_DIR)/summary.json: $(MODELS_DIR)/.done scripts/summarize.py | $(REPORT_DIR)
	$(RUN_AMR) python3 scripts/summarize.py \
		--models-dir $(MODELS_DIR) \
		--antibiotics $(ANTIBIOTICS) \
		--output $@

report:
	$(call need-samples,_report)

# Separate recipe lines, not prerequisites: analyze consumes what
# process-samples produces, and prerequisites of the same target are fair game
# for parallel execution under -j.
#
# The final cleanup removes everything under data/ except the assemblies, plus
# the feature matrices already folded into the trained models. See the
# disk-usage notes in README.md.
_report:
	$(MAKE) --no-print-directory process-samples
	$(MAKE) --no-print-directory analyze
	rm -rf $(READS_DIR) $(TRIMMED_DIR) $(DOWNSAMPLED_DIR) $(QC_DIR) \
	       $(KRAKEN_DIR) $(QUAST_DIR) $(AMR_DIR) $(CLEAN_DIR) \
	       $(FEATURES_DIR) $(SEQUENCES_DIR)
