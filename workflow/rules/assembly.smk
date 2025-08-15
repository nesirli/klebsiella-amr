rule spades_assembly:
    input:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz"
    output:
        "results/assemblies/{sample}/contigs.fasta"