process FASTP {
    tag "$run_accession"
    container 'quay.io/biocontainers/fastp:1.3.3--h43da1c4_0'
    cpus 4

    input:
    tuple val(run_accession), path(reads1), path(reads2)

    output:
    tuple val(run_accession), path("${run_accession}_trimmed_1.fastq.gz"), path("${run_accession}_trimmed_2.fastq.gz"), emit: trimmed
    tuple val(run_accession), path("${run_accession}_fastp.json"), emit: json
    tuple val(run_accession), path("${run_accession}_fastp.html"), emit: html

    script:
    """
    fastp \\
        --in1 ${reads1} \\
        --in2 ${reads2} \\
        --out1 ${run_accession}_trimmed_1.fastq.gz \\
        --out2 ${run_accession}_trimmed_2.fastq.gz \\
        --json ${run_accession}_fastp.json \\
        --html ${run_accession}_fastp.html \\
        --qualified_quality_phred 20 \\
        --length_required 50 \\
        --thread ${task.cpus}
    """
}
