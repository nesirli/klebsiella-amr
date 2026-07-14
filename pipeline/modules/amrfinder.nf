process AMRFINDER_DB {
    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.2.5--hf69ffd2_0'
    storeDir "${projectDir}/reference"

    output:
    path "amrfinderplus_db"

    script:
    """
    amrfinder_update -d amrfinderplus_db
    """
}

process AMRFINDER {
    tag "$run_accession"
    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.2.5--hf69ffd2_0'
    cpus 4

    input:
    tuple val(run_accession), path(assembly_fasta)
    path db_dir

    output:
    tuple val(run_accession), path("${run_accession}_amr.tsv"), emit: report

    script:
    """
    amrfinder \\
        --nucleotide ${assembly_fasta} \\
        --organism Klebsiella_pneumoniae \\
        --database ${db_dir}/latest \\
        --output ${run_accession}_amr.tsv \\
        --threads ${task.cpus}
    """
}
