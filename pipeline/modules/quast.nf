process QUAST {
    tag "$run_accession"
    container 'quay.io/biocontainers/quast:5.3.0--py39pl5321heaaa4ec_0'
    cpus 4

    input:
    tuple val(run_accession), path(assembly_fasta)

    output:
    tuple val(run_accession), path("${run_accession}_quast_report.tsv"), emit: report

    script:
    """
    quast.py \\
        --output-dir quast_out \\
        --threads ${task.cpus} \\
        ${assembly_fasta}

    mv quast_out/report.tsv ${run_accession}_quast_report.tsv
    """
}
