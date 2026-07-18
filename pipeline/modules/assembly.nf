process ASSEMBLY {
    tag "$run_accession"
    container 'quay.io/biocontainers/spades:4.2.0--h8d6e82b_1'
    cpus 6
    memory '10 GB'

    input:
    tuple val(run_accession), path(reads1), path(reads2)

    output:
    tuple val(run_accession), path("${run_accession}_assembled.fasta"), emit: fasta

    script:
    """
    spades.py \\
        --pe1-1 ${reads1} \\
        --pe1-2 ${reads2} \\
        --isolate \\
        --only-assembler \\
        -k 21,33,55 \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toGiga()} \\
        -o assembly_tmp

    mv assembly_tmp/contigs.fasta ${run_accession}_assembled.fasta
    rm -rf assembly_tmp
    """
}
