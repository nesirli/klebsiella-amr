nextflow.enable.dsl = 2

include { PARSE_METADATA } from './modules/metadata.nf'

workflow {
    PARSE_METADATA(file(params.metadata))
}
