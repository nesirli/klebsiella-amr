nextflow.enable.dsl = 2

include { PARSE_METADATA } from './modules/metadata.nf'
include { DOWNLOAD_READS } from './modules/download.nf'
include { FASTP }          from './modules/fastp.nf'

workflow {
    PARSE_METADATA(file(params.metadata))

    samples_ch = PARSE_METADATA.out[0]
        .mix(PARSE_METADATA.out[1])
        .splitCsv(header: true)
        .map { row -> row.run }
        .take(params.max_samples.toInteger())

    DOWNLOAD_READS(samples_ch)
    FASTP(DOWNLOAD_READS.out)
}
