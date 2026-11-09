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
// include { samplesheetToList         } from 'plugin/nf-schema'
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
    samples_ch = channel
        .fromPath(input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.sample
            ]

            tuple(
                meta,
                file(row.bam, checkIfExists: true),
                file(row.assembly, checkIfExists: true)
            )
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
        log.info "ONT methylation workflow initialised successfully"
    }
}