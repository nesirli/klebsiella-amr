# Auto-generated config export ------------------------------------------------
results/metadata/config.mk: config.yaml scripts/export_config.py
	mkdir -p results/metadata
	$(RUN_AMR) python3 scripts/export_config.py --config config.yaml --output $@
