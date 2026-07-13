process PARSE_METADATA {
    publishDir "${params.outdir}/metadata", mode: 'copy'

    input:
    path metadata

    output:
    path "train.csv"
    path "test.csv"

    script:
    def antibiotics = params.antibiotics.join(' ')
    def test_years  = params.test_years.join(' ')
    """
    uv run --project ${projectDir} python ${projectDir}/scripts/metadata.py \\
        --input ${metadata} \\
        --train-output train.csv \\
        --test-output test.csv \\
        --train-cutoff ${params.train_cutoff} \\
        --test-years ${test_years} \\
        --antibiotics ${antibiotics}
    """
}
