# Hyperparameter tuning (Optuna) ----------------------------------------------
# Each tuner writes {"hyperparameters": {...}, "tuning": {...}}; the matching
# trainer reads the first block via --params-input.
tune: $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/xgboost/tune_$(abx).json) \
      $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/lightgbm/tune_$(abx).json) \
      $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/nn/tune_$(abx).json)

$(MODELS_DIR)/xgboost/tune_%.json: $(FEATURES_DIR)/train_features.csv scripts/tune_xgboost.py scripts/common.py | $(MODELS_DIR)/xgboost
	$(RUN_AMR) python3 scripts/tune_xgboost.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--n-trials $(TUNE_TRIALS_XGB) \
		--n-splits $(TUNE_SPLITS) \
		--output $@

$(MODELS_DIR)/lightgbm/tune_%.json: $(FEATURES_DIR)/train_features.csv scripts/tune_lightgbm.py scripts/common.py | $(MODELS_DIR)/lightgbm
	$(RUN_AMR) python3 scripts/tune_lightgbm.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--n-trials $(TUNE_TRIALS_LGB) \
		--n-splits $(TUNE_SPLITS) \
		--output $@

$(MODELS_DIR)/nn/tune_%.json: $(FEATURES_DIR)/train_features.csv scripts/tune_nn.py scripts/mlp.py scripts/common.py | $(MODELS_DIR)/nn
	$(RUN_AMR) python3 scripts/tune_nn.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--n-trials $(TUNE_TRIALS_NN) \
		--n-splits $(TUNE_SPLITS) \
		--output $@

# Models ----------------------------------------------------------------------
models:
	@if [ -z "$(SAMPLES)" ]; then \
		$(MAKE) --no-print-directory metadata MAX_SAMPLES=$(MAX_SAMPLES) && $(MAKE) --no-print-directory models MAX_SAMPLES=$(MAX_SAMPLES); \
	else \
		$(MAKE) --no-print-directory $(MODELS_DIR)/.done; \
	fi

$(MODELS_DIR)/.done: $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/xgboost/$(abx)_metrics.json) \
                     $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/lightgbm/$(abx)_metrics.json) \
                     $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/nn/$(abx)_metrics.json)
	@touch $@

$(MODELS_DIR)/xgboost/%_metrics.json: $(FEATURES_DIR)/train_features.csv $(FEATURES_DIR)/test_features.csv scripts/train_xgboost.py scripts/common.py $(MODELS_DIR)/xgboost/tune_%.json | $(MODELS_DIR)/xgboost
	$(RUN_AMR) python3 scripts/train_xgboost.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--test-features $(FEATURES_DIR)/test_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--params-input $(MODELS_DIR)/xgboost/tune_$*.json \
		--model-output $(MODELS_DIR)/xgboost/$*_model.json \
		--params-output $(MODELS_DIR)/xgboost/$*_params.json \
		--metrics-output $@ \
		--predictions-output $(MODELS_DIR)/xgboost/$*_predictions.csv \
		--importance-output $(MODELS_DIR)/xgboost/$*_importance.csv \
		--importance-plot-output $(MODELS_DIR)/xgboost/$*_shap.png

$(MODELS_DIR)/lightgbm/%_metrics.json: $(FEATURES_DIR)/train_features.csv $(FEATURES_DIR)/test_features.csv scripts/train_lightgbm.py scripts/common.py $(MODELS_DIR)/lightgbm/tune_%.json | $(MODELS_DIR)/lightgbm
	$(RUN_AMR) python3 scripts/train_lightgbm.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--test-features $(FEATURES_DIR)/test_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--params-input $(MODELS_DIR)/lightgbm/tune_$*.json \
		--model-output $(MODELS_DIR)/lightgbm/$*_model.txt \
		--params-output $(MODELS_DIR)/lightgbm/$*_params.json \
		--metrics-output $@ \
		--predictions-output $(MODELS_DIR)/lightgbm/$*_predictions.csv \
		--importance-output $(MODELS_DIR)/lightgbm/$*_importance.csv \
		--importance-plot-output $(MODELS_DIR)/lightgbm/$*_shap.png

$(MODELS_DIR)/nn/%_metrics.json: $(FEATURES_DIR)/train_features.csv $(FEATURES_DIR)/test_features.csv scripts/train_nn.py scripts/mlp.py scripts/common.py $(MODELS_DIR)/nn/tune_%.json | $(MODELS_DIR)/nn
	$(RUN_AMR) python3 scripts/train_nn.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--test-features $(FEATURES_DIR)/test_features.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--params-input $(MODELS_DIR)/nn/tune_$*.json \
		--model-output $(MODELS_DIR)/nn/$*_model.pt \
		--params-output $(MODELS_DIR)/nn/$*_params.json \
		--metrics-output $@ \
		--predictions-output $(MODELS_DIR)/nn/$*_predictions.csv \
		--importance-output $(MODELS_DIR)/nn/$*_importance.csv \
		--importance-plot-output $(MODELS_DIR)/nn/$*_importance.png

# DNABERT-2 (optional) --------------------------------------------------------
dnabert:
	@if [ -z "$(SAMPLES)" ]; then \
		$(MAKE) --no-print-directory metadata MAX_SAMPLES=$(MAX_SAMPLES) && $(MAKE) --no-print-directory dnabert MAX_SAMPLES=$(MAX_SAMPLES); \
	else \
		$(MAKE) --no-print-directory $(MODELS_DIR)/dnabert/.done; \
	fi

$(MODELS_DIR)/dnabert/.done: $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/dnabert/$(abx)_metrics.json) | $(MODELS_DIR)/dnabert
	@touch $@

$(MODELS_DIR)/dnabert/%_metrics.json: $(SEQUENCES_DIR)/train_sequences.csv $(SEQUENCES_DIR)/test_sequences.csv scripts/train_dnabert.py scripts/common.py | $(MODELS_DIR)/dnabert
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
