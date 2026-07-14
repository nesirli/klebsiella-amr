process DOWNSAMPLE {
    tag "$run_accession"
    container 'quay.io/biocontainers/seqtk:1.5--h577a1d6_1'
    cpus 2

    input:
    tuple val(run_accession), path(reads1), path(reads2), path(fastp_json)

    output:
    tuple val(run_accession), path("${run_accession}_downsampled_1.fastq.gz"), path("${run_accession}_downsampled_2.fastq.gz"), emit: reads

    script:
    """
    total_bases=\$(grep -A2 '"after_filtering"' ${fastp_json} | grep total_bases | grep -oE '[0-9]+')
    fraction=\$(awk -v tb="\$total_bases" -v gs="${params.genome_size}" -v tc="${params.target_coverage}" \\
        'BEGIN { cov = tb / gs; f = tc / cov; if (f > 1) f = 1; printf "%.6f", f }')

    seqtk sample -s42 ${reads1} "\$fraction" | gzip > ${run_accession}_downsampled_1.fastq.gz
    seqtk sample -s42 ${reads2} "\$fraction" | gzip > ${run_accession}_downsampled_2.fastq.gz
    """
}
