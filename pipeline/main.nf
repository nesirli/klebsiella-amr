nextflow.enable.dsl = 2

include { PARSE_METADATA }        from './modules/metadata.nf'
include { DOWNLOAD_READS }        from './modules/download.nf'
include { FASTP }                 from './modules/fastp.nf'
include { KRAKEN2_DB; KRAKEN2 }   from './modules/kraken2.nf'
include { DOWNSAMPLE }            from './modules/downsample.nf'

workflow {
    main:
    PARSE_METADATA(file(params.metadata))

    samples_ch = PARSE_METADATA.out.train
        .mix(PARSE_METADATA.out.test)
        .splitCsv(header: true)
        .map { row -> row.run }
        .take(params.max_samples.toInteger())

    DOWNLOAD_READS(samples_ch)
    FASTP(DOWNLOAD_READS.out)

    KRAKEN2_DB()
    KRAKEN2(FASTP.out.trimmed, KRAKEN2_DB.out.first())

    DOWNSAMPLE(FASTP.out.trimmed.join(FASTP.out.json))

    publish:
    metadata_train = PARSE_METADATA.out.train
    metadata_test  = PARSE_METADATA.out.test
    qc_json        = FASTP.out.json
    qc_html        = FASTP.out.html
    kraken2_report = KRAKEN2.out.report
    downsampled    = DOWNSAMPLE.out.reads
}

output {
    metadata_train {
        path 'metadata'
    }
    metadata_test {
        path 'metadata'
    }
    qc_json {
        path 'qc'
    }
    qc_html {
        path 'qc'
    }
    kraken2_report {
        path 'kraken2'
    }
    downsampled {
        path 'downsampled'
    }
}
