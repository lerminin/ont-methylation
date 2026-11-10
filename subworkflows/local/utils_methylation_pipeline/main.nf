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

            def meta = [id: row.sample]

            def fasta = file(row.fasta, checkIfExists: true)
            def bam   = file(row.bam, checkIfExists: true)

            def bins = params.meta
                ? file(row.bins_folder, checkIfExists: true)
                : null

            tuple(meta, fasta, bam, bins)
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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP MESSAGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// def helpMSG() {
//         c_green = "\033[0;32m";
//         c_reset = "\033[0m";
//         c_yellow = "\033[0;33m";
//         c_blue = "\033[0;34m";
//         c_red = "\033[0;31m";
//         c_dim = "\033[2m";
//         log.info """
//         ____________________________________________________________________________________________

//         Nextflow Pipeline for Methylated Motif Extraction and Statistical Analysis from ONT data.

//         ${c_yellow}Usage example:${c_reset}
//         nextflow run valegale/ONT_methylation -r 0.0.1 --fasta '*.fasta' --bam '*.bam' 

//         Use the following commands to check for latest pipeline versions:
        
//         nextflow pull valegale/ONT_methylation
//         nextflow info valegale/ONT_methylation

//         ${c_yellow}Input${c_reset}
//         ${c_green} --fasta ${c_reset}           '*.fasta'       -> one genome/assembly per file
//         ${c_green} --bam ${c_reset}             '*.bam'         -> one sorted BAM matching one reference FASTA

//         ${c_dim}  change above input to csv:${c_reset} ${c_green}--list ${c_reset}

//         ${c_yellow}IMPORTANT:${c_reset} Unless ${c_green}--list${c_reset} is used, the ${c_yellow}basename${c_reset} of the FASTA and BAM files must match
//         (e.g., sample1.fasta <-> sample1.bam). 

//         ${c_yellow}General Options:${c_reset}
//         --cores             Max cores per process for local use [default: $params.cores]
//         --max_cores         Max cores (in total) for local use [default: $params.max_cores]
//         --memory            Max memory for local use [default: $params.memory]
//         --outdir            Name of the result folder [default: $params.outdir]

//         ${c_yellow}Additional Options:${c_reset}${c_reset}
//         --filter_threshold_modkit             Filter threshold for modkit [default: $params.filter_threshold_modkit]
//         --automatic_threshold_modkit          Enable automatic estimation of the filter threshold by modkit.
//                                             When true, modkit will determine an optimal threshold from the data and the value of --filter_threshold_modkit will be ignored.
//                                             [default: $params.automatic_threshold_modkit]
//         --percent_cutoff_modification_table   Minimum methylation percentage required for genome positions to be reported in the modification tables [default: $params.percent_cutoff_modification_table].
    
//         ${c_dim}Nextflow options:
//         -with-report rep.html    cpu / ram usage (may cause errors)
//         -with-dag chart.html     generates a flowchart for the process tree
//         -with-timeline time.html timeline (may cause errors)
//         -resume                  resume a previous calculation w/o recalculating everything (needs the same run command and work dir!)

//         ${c_yellow}Caching:${c_reset}
//         --singularityCacheDir   Location for storing the Singularity images [default: $params.singularityCacheDir]
//         -w                      Working directory for all intermediate results [default: work] 

//         ${c_yellow}Execution/Engine profiles:${c_reset}
//         The pipeline supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}docker${c_reset}
        
//         ${c_green}Executer${c_reset} (choose one):
//         local
//         slurm
        
//         ${c_blue}Engines${c_reset} (choose one):
//         docker
//         singularity
        
//         Per default: -profile local,docker is executed (-profile standard).
        
//         ${c_reset}
//         """.stripIndent()
//     }