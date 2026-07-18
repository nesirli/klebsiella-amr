process BUILD_FEATURES {
    container 'ghcr.io/astral-sh/uv:python3.12-bookworm-slim'

    input:
    path amr_files
    path train_metadata
    path test_metadata

    output:
    path "train_features.csv", emit: train
    path "test_features.csv",  emit: test

    script:
    """
    uv run "\$(command -v build_features.py)" \\
        --amr-files ${amr_files} \\
        --train-metadata ${train_metadata} \\
        --test-metadata ${test_metadata} \\
        --train-output train_features.csv \\
        --test-output test_features.csv
    """
}
