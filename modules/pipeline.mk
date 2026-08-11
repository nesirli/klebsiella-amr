# Active sample tracking ------------------------------------------------------
SAMPLES_ACTIVE       = $(sort $(patsubst $(QC_DIR)/%_fastp.json,%,$(wildcard $(QC_DIR)/*_fastp.json)))
SAMPLES_WITH_READS   = $(sort $(patsubst $(READS_DIR)/%_1.fastq.gz,%,$(wildcard $(READS_DIR)/*_1.fastq.gz)))

# Reads -----------------------------------------------------------------------
$(READS_DIR)/.done: $(patsubst %,$(READS_DIR)/%_1.fastq.gz,$(SAMPLES_ACTIVE)) | $(READS_DIR)
	@touch $@

$(READS_DIR)/%_1.fastq.gz $(READS_DIR)/%_2.fastq.gz &: scripts/download_reads.py | $(READS_DIR)
	$(RUN_AMR) python3 scripts/download_reads.py \
		--accession $* \
		--out1 $(READS_DIR)/$*_1.fastq.gz \
		--out2 $(READS_DIR)/$*_2.fastq.gz \
		--allow-skip

# Trimming / fastp ------------------------------------------------------------
$(TRIMMED_DIR)/.done: $(patsubst %,$(QC_DIR)/%_fastp.json,$(SAMPLES_ACTIVE)) | $(TRIMMED_DIR)
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

$(QC_DIR)/.done: $(patsubst %,$(QC_DIR)/%_fastp.json,$(SAMPLES_ACTIVE)) | $(QC_DIR)
	@touch $@

# Kraken2 ---------------------------------------------------------------------
$(KRAKEN_DB)/taxo.k2d:
	mkdir -p $(KRAKEN_DB)
	curl -L -o $(KRAKEN_DB)/k2_standard_08gb_20250402.tar.gz \
		https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz
	tar -xzf $(KRAKEN_DB)/k2_standard_08gb_20250402.tar.gz -C $(KRAKEN_DB)
	rm $(KRAKEN_DB)/k2_standard_08gb_20250402.tar.gz

$(KRAKEN_DIR)/.done: $(patsubst %,$(KRAKEN_DIR)/%_report.txt,$(SAMPLES_ACTIVE)) | $(KRAKEN_DIR)
	@touch $@

$(KRAKEN_DIR)/%_report.txt: $(TRIMMED_DIR)/%_1.fastq.gz $(TRIMMED_DIR)/%_2.fastq.gz $(KRAKEN_DB)/taxo.k2d | $(KRAKEN_DIR)
	$(RUN_BIOINFO) kraken2 --db $(KRAKEN_DB) \
		--paired --gzip-compressed --memory-mapping \
		--report $@ --output /dev/null \
		$(TRIMMED_DIR)/$*_1.fastq.gz $(TRIMMED_DIR)/$*_2.fastq.gz

# QUAST -----------------------------------------------------------------------
$(QUAST_DIR)/.done: $(patsubst %,$(QUAST_DIR)/%/report.tsv,$(SAMPLES_ACTIVE)) | $(QUAST_DIR)
	@touch $@

$(QUAST_DIR)/%/report.tsv: $(ASSEMBLY_DIR)/%_assembled.fasta | $(QUAST_DIR)
	$(RUN_BIOINFO) quast.py \
		--output-dir $(QUAST_DIR)/$* \
		--threads $(QUAST_THREADS) \
		$<

# Assembly --------------------------------------------------------------------
$(DOWNSAMPLED_DIR)/.done: $(patsubst %,$(DOWNSAMPLED_DIR)/%_1.fastq.gz,$(SAMPLES_ACTIVE)) | $(DOWNSAMPLED_DIR)
	@touch $@

$(DOWNSAMPLED_DIR)/%_1.fastq.gz $(DOWNSAMPLED_DIR)/%_2.fastq.gz &: $(TRIMMED_DIR)/%_1.fastq.gz $(TRIMMED_DIR)/%_2.fastq.gz $(QC_DIR)/%_fastp.json | $(DOWNSAMPLED_DIR)
	fraction=$$($(RUN_AMR) python3 -c "import json; d=json.load(open('$(QC_DIR)/$*_fastp.json')); cov=d['summary']['after_filtering']['total_bases']/$(GENOME_SIZE); print(min(0.999999, $(TARGET_COVERAGE)/cov))"); \
	$(RUN_BIOINFO) seqtk sample -s42 $(TRIMMED_DIR)/$*_1.fastq.gz $$fraction | gzip > $(DOWNSAMPLED_DIR)/$*_1.fastq.gz; \
	$(RUN_BIOINFO) seqtk sample -s42 $(TRIMMED_DIR)/$*_2.fastq.gz $$fraction | gzip > $(DOWNSAMPLED_DIR)/$*_2.fastq.gz

$(ASSEMBLY_DIR)/.done: $(patsubst %,$(ASSEMBLY_DIR)/%_assembled.fasta,$(SAMPLES_ACTIVE)) | $(ASSEMBLY_DIR)
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

# AMRFinderPlus ---------------------------------------------------------------
$(AMR_DIR)/.done: $(patsubst %,$(AMR_DIR)/%_amr.tsv,$(SAMPLES_ACTIVE)) | $(AMR_DIR)
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

# Per-sample cleanup ----------------------------------------------------------
$(CLEAN_DIR)/%_cleaned: $(AMR_DIR)/%_amr.tsv $(KRAKEN_DIR)/%_report.txt $(QUAST_DIR)/%/report.tsv | $(CLEAN_DIR)
	@rm -f $(READS_DIR)/$*_1.fastq.gz $(READS_DIR)/$*_2.fastq.gz \
	      $(TRIMMED_DIR)/$*_1.fastq.gz $(TRIMMED_DIR)/$*_2.fastq.gz \
	      $(DOWNSAMPLED_DIR)/$*_1.fastq.gz $(DOWNSAMPLED_DIR)/$*_2.fastq.gz
	@touch $@
