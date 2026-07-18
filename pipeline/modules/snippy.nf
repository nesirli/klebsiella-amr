process SNIPPY {
    tag "$run_accession"
    container 'quay.io/biocontainers/snippy:4.6.0--hdfd78af_6'
    cpus 4
    memory '4 GB'

    input:
    tuple val(run_accession), path(reads1), path(reads2)
    path reference_fasta

    output:
    tuple val(run_accession), path("${run_accession}/snps.vcf"), emit: vcf

    script:
    """
    snippy \\
        --outdir ${run_accession} \\
        --reference ${reference_fasta} \\
        --R1 ${reads1} \\
        --R2 ${reads2} \\
        --cpus ${task.cpus} \\
        --cleanup \\
        --force
    """
}
