process DOWNLOAD_READS {
    tag "$run_accession"
    maxForks 4

    input:
    val run_accession

    output:
    tuple val(run_accession), path("${run_accession}_1.fastq.gz"), path("${run_accession}_2.fastq.gz")

    script:
    """
    urls=\$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run_accession}&result=read_run&fields=fastq_ftp&format=tsv" | tail -n +2 | cut -f2)
    url1=\$(echo "\$urls" | cut -d';' -f1)
    url2=\$(echo "\$urls" | cut -d';' -f2)
    curl --retry 5 --retry-delay 5 --retry-all-errors -C - -L -o ${run_accession}_1.fastq.gz "https://\$url1"
    curl --retry 5 --retry-delay 5 --retry-all-errors -C - -L -o ${run_accession}_2.fastq.gz "https://\$url2"
    """
}
