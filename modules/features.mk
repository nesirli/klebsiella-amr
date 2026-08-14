# Features (tabular) ----------------------------------------------------------
# Both splits come out of one invocation, so they are a grouped target ('&:'):
# make runs the recipe once and knows it produced all of them.
$(FEATURES_DIR)/train_features.csv $(FEATURES_DIR)/test_features.csv &: \
		$(patsubst %,$(CLEAN_DIR)/%_cleaned,$(SAMPLES_ACTIVE)) \
		results/metadata/train.csv results/metadata/test.csv \
		scripts/build_features.py | $(FEATURES_DIR)
	$(RUN_AMR) python3 scripts/build_features.py \
		--amr-files $(AMR_DIR)/*_amr.tsv \
		--train-metadata results/metadata/train.csv \
		--test-metadata results/metadata/test.csv \
		--train-output $(FEATURES_DIR)/train_features.csv \
		--test-output $(FEATURES_DIR)/test_features.csv

# Sequences (for optional DNABERT-2) ------------------------------------------
$(SEQUENCES_DIR)/train_sequences.csv $(SEQUENCES_DIR)/test_sequences.csv &: \
		$(patsubst %,$(CLEAN_DIR)/%_cleaned,$(SAMPLES_ACTIVE)) \
		results/metadata/train.csv results/metadata/test.csv \
		scripts/build_sequences.py | $(SEQUENCES_DIR)
	$(RUN_AMR) python3 scripts/build_sequences.py \
		--seq-files $(AMR_DIR)/*_amr_genes.fna \
		--train-metadata results/metadata/train.csv \
		--test-metadata results/metadata/test.csv \
		--train-output $(SEQUENCES_DIR)/train_sequences.csv \
		--test-output $(SEQUENCES_DIR)/test_sequences.csv
