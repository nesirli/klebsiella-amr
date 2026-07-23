nextflow.enable.dsl = 2

include { PARSE_METADATA }        from './modules/metadata.nf'
include { DOWNLOAD_READS }        from './modules/download.nf'
include { FASTP }                 from './modules/fastp.nf'
include { KRAKEN2_DB; KRAKEN2 }   from './modules/kraken2.nf'
include { DOWNSAMPLE }            from './modules/downsample.nf'
include { ASSEMBLY }              from './modules/assembly.nf'
include { QUAST }                 from './modules/quast.nf'
include { MULTIQC as MULTIQC_FASTP }   from './modules/multiqc.nf'
include { MULTIQC as MULTIQC_KRAKEN2 } from './modules/multiqc.nf'
include { MULTIQC as MULTIQC_QUAST }   from './modules/multiqc.nf'
include { AMRFINDER_DB; AMRFINDER }    from './modules/amrfinder.nf'
include { REFERENCE_GENOME }           from './modules/reference.nf'
include { SNIPPY }                     from './modules/snippy.nf'
include { BUILD_FEATURES }             from './modules/features.nf'
include { TRAIN_MODEL }                from './modules/train_model.nf'

workflow {
    main:
    PARSE_METADATA(file(params.metadata))

    // Cap train and test cohorts independently so a small dev run always
    // spans both year-splits (needed for the model stage to have training
    // data) rather than clustering into one via mix()+take(). Selection is
    // deterministic per cohort since splitCsv preserves metadata row order.
    def max_samples = params.max_samples.toInteger()

    train_ids = PARSE_METADATA.out.train
        .splitCsv(header: true)
        .map { row -> row.run }
        .take(max_samples)

    test_ids = PARSE_METADATA.out.test
        .splitCsv(header: true)
        .map { row -> row.run }
        .take(max_samples)

    samples_ch = train_ids.mix(test_ids)

    DOWNLOAD_READS(samples_ch)
    FASTP(DOWNLOAD_READS.out)

    KRAKEN2_DB()
    KRAKEN2(FASTP.out.trimmed, KRAKEN2_DB.out.first())

    DOWNSAMPLE(FASTP.out.trimmed.join(FASTP.out.json))
    ASSEMBLY(DOWNSAMPLE.out.reads)
    QUAST(ASSEMBLY.out.fasta)

    MULTIQC_FASTP('fastp',     FASTP.out.json.map { it[1] }.collect())
    MULTIQC_KRAKEN2('kraken2', KRAKEN2.out.report.map { it[1] }.collect())
    MULTIQC_QUAST('quast',     QUAST.out.report.map { it[1] }.collect())

    AMRFINDER_DB()
    AMRFINDER(ASSEMBLY.out.fasta, AMRFINDER_DB.out.first())

    REFERENCE_GENOME()
    SNIPPY(FASTP.out.trimmed, REFERENCE_GENOME.out.first())

    BUILD_FEATURES(
        AMRFINDER.out.report.map { it[1] }.collect(),
        PARSE_METADATA.out.train,
        PARSE_METADATA.out.test
    )

    antibiotics_ch = Channel.fromList(params.antibiotics)
    TRAIN_MODEL(antibiotics_ch, BUILD_FEATURES.out.train, BUILD_FEATURES.out.test)

    publish:
    metadata_train   = PARSE_METADATA.out.train
    metadata_test    = PARSE_METADATA.out.test
    qc_json          = FASTP.out.json
    qc_html          = FASTP.out.html
    kraken2_report   = KRAKEN2.out.report
    downsampled      = DOWNSAMPLE.out.reads
    assembly         = ASSEMBLY.out.fasta
    quast_report     = QUAST.out.report
    fastp_multiqc    = MULTIQC_FASTP.out.report
    kraken2_multiqc  = MULTIQC_KRAKEN2.out.report
    quast_multiqc    = MULTIQC_QUAST.out.report
    amr_report        = AMRFINDER.out.report
    snippy_vcf        = SNIPPY.out.vcf
    train_features    = BUILD_FEATURES.out.train
    test_features     = BUILD_FEATURES.out.test
    model_metrics     = TRAIN_MODEL.out.metrics
    model_predictions = TRAIN_MODEL.out.predictions
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
    assembly {
        path 'assembly'
    }
    quast_report {
        path 'quast'
    }
    fastp_multiqc {
        path 'multiqc'
    }
    kraken2_multiqc {
        path 'multiqc'
    }
    quast_multiqc {
        path 'multiqc'
    }
    amr_report {
        path 'amr'
    }
    snippy_vcf {
        path 'snippy'
    }
    train_features {
        path 'features'
    }
    test_features {
        path 'features'
    }
    model_metrics {
        path 'models'
    }
    model_predictions {
        path 'models'
    }
}
