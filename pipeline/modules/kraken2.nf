process KRAKEN2_DB {
    storeDir "${projectDir}/reference"

    output:
    path "kraken2_db"

    script:
    """
    mkdir -p kraken2_db
    curl -L -o k2.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz
    tar -xzf k2.tar.gz -C kraken2_db
    rm k2.tar.gz
    """
}

process KRAKEN2 {
    tag "$run_accession"
    container 'quay.io/biocontainers/kraken2:2.1.5--pl5321h077b44d_0'
    cpus 4

    input:
    tuple val(run_accession), path(reads1), path(reads2)
    path db_dir

    output:
    tuple val(run_accession), path("${run_accession}_kraken2_report.txt"), emit: report

    script:
    """
    kraken2 \\
        --db ${db_dir} \\
        --paired \\
        --gzip-compressed \\
        --report ${run_accession}_kraken2_report.txt \\
        --output /dev/null \\
        --threads ${task.cpus} \\
        ${reads1} ${reads2}
    """
}
