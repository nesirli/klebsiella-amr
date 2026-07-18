process MULTIQC {
    tag "$report_name"
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'

    input:
    val report_name
    path qc_files, stageAs: 'input*/*'

    output:
    path "${report_name}_multiqc_report.html", emit: report

    script:
    """
    multiqc . --filename ${report_name}_multiqc_report.html --force
    """
}
