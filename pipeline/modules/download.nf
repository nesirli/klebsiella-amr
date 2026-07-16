process DOWNLOAD_READS {
    tag "$run_accession"
    container 'ghcr.io/astral-sh/uv:python3.12-bookworm-slim'
    maxForks 4

    input:
    val run_accession

    output:
    tuple val(run_accession), path("${run_accession}_1.fastq.gz"), path("${run_accession}_2.fastq.gz")

    script:
    """
    python3 "\$(command -v download_reads.py)" \\
        --accession ${run_accession} \\
        --out1 ${run_accession}_1.fastq.gz \\
        --out2 ${run_accession}_2.fastq.gz
    """
}
