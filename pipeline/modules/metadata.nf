process PARSE_METADATA {
    container 'ghcr.io/astral-sh/uv:python3.12-bookworm-slim'

    input:
    path metadata

    output:
    path "train.csv", emit: train
    path "test.csv",  emit: test

    script:
    def antibiotics = params.antibiotics.join(' ')
    def test_years  = params.test_years.join(' ')
    """
    uv run "\$(command -v metadata.py)" \\
        --input ${metadata} \\
        --train-output train.csv \\
        --test-output test.csv \\
        --train-cutoff ${params.train_cutoff} \\
        --test-years ${test_years} \\
        --antibiotics ${antibiotics}
    """
}
