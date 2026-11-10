//
// Subworkflow with functionality specific to the lerminin/ont-methylation pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
// include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
// include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
// include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
// include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
// include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    input

    main:

    if (!input) {
        error "Please provide a sample sheet using --input"
    }

    def required_columns = params.meta
    ? ['sample', 'bam', 'assembly', 'bins_folder']
    : ['sample', 'bam', 'assembly']

    samples_ch = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->

            required_columns.each { col ->
                if (!row[col]?.toString()?.trim()) {
                    error "Missing required value for '${col}' in samplesheet"
                }
            }

            def meta = [id: row.sample]

            def bam   = file(row.bam, checkIfExists: true)
            def assembly = file(row.assembly, checkIfExists: true)

            def bins = params.meta
                ? file(row.bins_folder, checkIfExists: true)
                : null

            tuple(meta, bam, assembly, bins)
        }

    emit:
    samplesheet = samples_ch
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    main:
    workflow.onComplete {
        //log.info "Pipeline completed: ${workflow.success}"
        //log.info "Duration: ${workflow.duration}"
        log.info "ONT methylation workflow completed successfully"
    }
}
