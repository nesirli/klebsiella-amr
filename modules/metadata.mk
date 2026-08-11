# Metadata --------------------------------------------------------------------
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
