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
	$(call need-samples,$(MODELS_DIR)/.done)

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
	$(call need-samples,$(MODELS_DIR)/dnabert/.done)

$(MODELS_DIR)/dnabert/.done: $(foreach abx,$(ANTIBIOTICS),$(MODELS_DIR)/dnabert/$(abx)_metrics.json) | $(MODELS_DIR)/dnabert
	@touch $@

$(MODELS_DIR)/dnabert/tune_%.json: $(SEQUENCES_DIR)/train_sequences.csv scripts/tune_dnabert.py scripts/train_dnabert.py scripts/common.py | $(MODELS_DIR)/dnabert
	$(RUN_AMR) python3 scripts/tune_dnabert.py \
		--train-features $(SEQUENCES_DIR)/train_sequences.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--n-trials $(TUNE_TRIALS_DNABERT) \
		--n-splits $(TUNE_SPLITS_DNABERT) \
		--max-train $(TUNE_MAX_TRAIN_DNABERT) \
		--output $@

$(MODELS_DIR)/dnabert/%_metrics.json: $(SEQUENCES_DIR)/train_sequences.csv $(SEQUENCES_DIR)/test_sequences.csv scripts/train_dnabert.py scripts/common.py $(MODELS_DIR)/dnabert/tune_%.json | $(MODELS_DIR)/dnabert
	$(RUN_AMR) python3 scripts/train_dnabert.py \
		--train-features $(SEQUENCES_DIR)/train_sequences.csv \
		--test-features $(SEQUENCES_DIR)/test_sequences.csv \
		--antibiotic $* \
		--all-antibiotics $(ANTIBIOTICS) \
		--params-input $(MODELS_DIR)/dnabert/tune_$*.json \
		--model-output $(MODELS_DIR)/dnabert/$*_model.pt \
		--params-output $(MODELS_DIR)/dnabert/$*_params.json \
		--metrics-output $@ \
		--predictions-output $(MODELS_DIR)/dnabert/$*_predictions.csv \
		--importance-output $(MODELS_DIR)/dnabert/$*_importance.csv \
		--importance-plot-output $(MODELS_DIR)/dnabert/$*_importance.png

# Serving -----------------------------------------------------------------------
# The manifest freezes the feature layout and thresholds the app needs. It must
# be built while the feature matrices still exist: 'make report' deletes them.
manifest: $(MODELS_DIR)/manifest.json

$(MODELS_DIR)/manifest.json: $(FEATURES_DIR)/train_features.csv $(MODELS_DIR)/.done scripts/export_manifest.py
	$(RUN_AMR) python3 scripts/export_manifest.py \
		--train-features $(FEATURES_DIR)/train_features.csv \
		--models-dir $(MODELS_DIR) \
		--antibiotics $(ANTIBIOTICS) \
		--amrfinder-db $(AMRFINDER_DB) \
		--output $@

# Copy just what the container serves into app/artifacts: the manifest plus the
# two tree models. Small enough to commit (~2 MB), which keeps `docker build
# app/` self-contained -- Coolify needs no access to results/.
app-artifacts: $(MODELS_DIR)/manifest.json
	@mkdir -p app/artifacts/xgboost app/artifacts/lightgbm
	cp $(MODELS_DIR)/manifest.json app/artifacts/manifest.json
	cp $(MODELS_DIR)/xgboost/*_model.json app/artifacts/xgboost/
	cp $(MODELS_DIR)/lightgbm/*_model.txt app/artifacts/lightgbm/
	@echo 'app/artifacts ready'

app: app-artifacts
	$(RUN_AMR) streamlit run app/app.py
