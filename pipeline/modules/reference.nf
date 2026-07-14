process REFERENCE_GENOME {
    storeDir "${projectDir}/reference"

    output:
    path "reference_genome.fasta"

    script:
    """
    curl -L -o reference_genome.fasta.gz ${params.reference_url}
    gunzip reference_genome.fasta.gz
    """
}
