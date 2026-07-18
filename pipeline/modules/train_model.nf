process TRAIN_MODEL {
    tag "$antibiotic"
    container 'ghcr.io/astral-sh/uv:python3.12-bookworm-slim'

    input:
    val antibiotic
    path train_features
    path test_features

    output:
    tuple val(antibiotic), path("${antibiotic}_metrics.json"), emit: metrics
    tuple val(antibiotic), path("${antibiotic}_predictions.csv"), emit: predictions

    script:
    def all_antibiotics = params.antibiotics.join(' ')
    """
    uv run "\$(command -v train_model.py)" \\
        --train-features ${train_features} \\
        --test-features ${test_features} \\
        --antibiotic ${antibiotic} \\
        --all-antibiotics ${all_antibiotics} \\
        --metrics-output ${antibiotic}_metrics.json \\
        --predictions-output ${antibiotic}_predictions.csv
    """
}
