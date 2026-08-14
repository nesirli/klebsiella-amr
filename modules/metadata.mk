# Metadata --------------------------------------------------------------------
METADATA_STAMP := results/metadata/.max-samples

metadata: results/metadata/.done

results/metadata/train.csv results/metadata/test.csv results/metadata/samples.mk results/metadata/samples.txt: results/metadata/.done
	@true

# MAX_SAMPLES is a command-line knob, not a file, so make cannot see it change.
# Recording it in a stamp file that is only rewritten when the value differs
# makes 'make metadata MAX_SAMPLES=20' re-split after a run at 5, while leaving
# the split alone when the value is unchanged.
$(METADATA_STAMP): force
	@mkdir -p $(@D)
	@echo '$(MAX_SAMPLES)' | cmp -s - $@ || echo '$(MAX_SAMPLES)' > $@

force:

results/metadata/.done: config.yaml scripts/metadata.py results/metadata/config.mk \
                        $(METADATA) $(METADATA_STAMP)
	mkdir -p results/metadata
	$(RUN_AMR) python3 scripts/metadata.py \
		--config config.yaml \
		--train-output results/metadata/train.csv \
		--test-output results/metadata/test.csv \
		--samples-output results/metadata/samples.txt \
		--max-samples $(MAX_SAMPLES)
	touch $@
