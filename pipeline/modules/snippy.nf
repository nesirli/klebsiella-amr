process SNIPPY {
    tag "$run_accession"
    container 'quay.io/biocontainers/snippy:4.6.0--hdfd78af_6'
    cpus 4
    memory '8 GB'

    input:
    tuple val(run_accession), path(reads1), path(reads2)
    path reference_fasta

    output:
    tuple val(run_accession), path("${run_accession}/snps.vcf"), emit: vcf

    script:
    // Nextflow enforces `memory` as a hard docker --memory cgroup limit.
    // Snippy's default sizes its two concurrent `samtools sort` buffers to
    // ~all available RAM, which overshoots any container cap and gets the
    // samtools OOM-killed mid-pipe ("error reading header"). Pin --ram so
    // snippy's internal budget stays comfortably under the declared limit.
    """
    snippy \\
        --outdir ${run_accession} \\
        --reference ${reference_fasta} \\
        --R1 ${reads1} \\
        --R2 ${reads2} \\
        --cpus ${task.cpus} \\
        --ram 6 \\
        --cleanup \\
        --force
    """
}
