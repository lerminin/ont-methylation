/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM2FASTQ; ZIPFASTQ; MINIMAP2; SPLIT_BAM_BY_BIN } from '../modules/map_index_bam.nf'
include { MODKIT_PILEUP; MODKIT_PILEUP_BEDGRAPHS; MODKIT_FIND_MOTIFS; CUSTOM_BEDGRAPHS; PUBLISH_RESULTS_META; PUBLISH_RESULTS_MOTIFS_META; PUBLISH_RESULTS; PUBLISH_RESULTS_MOTIFS} from '../modules/modkit.nf'
include { COMPUTE_STATISTICS } from '../modules/statistics.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INITIALIZE CHANNELS FROM PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ONT_METHYLATION {
    take:
    ch_samplesheet

    main:
    fastq_ch = BAM2FASTQ(ch_samplesheet)
    
    // ZIPFASTQ(fastq_ch)

    mapped_bams_ch = MINIMAP2(BAM2FASTQ.out)

    pileup_ch = MODKIT_PILEUP(mapped_bams_ch)
    
    pileup_bedgraphs_ch = MODKIT_PILEUP_BEDGRAPHS(mapped_bams_ch)
    
    custom_bedgraphs_ch = CUSTOM_BEDGRAPHS(pileup_ch)

    motifs_ch = MODKIT_FIND_MOTIFS(pileup_ch)

    statistics_ch = COMPUTE_STATISTICS(pileup_ch)

    publish_input = pileup_ch.join(
            pileup_bedgraphs_ch.map { meta, bedgraphs, assembly ->
                tuple(meta, bedgraphs)
            }
        )
        .join(custom_bedgraphs_ch)
        .join(statistics_ch)

    publish_motifs_input = pileup_ch.join(motifs_ch)

    PUBLISH_RESULTS(publish_input)
    PUBLISH_RESULTS_MOTIFS(publish_motifs_input)
}

  
// }

// workflow HELP {

//     // terminal prints
//     println " "
//     println "\u001B[32mProfile: $workflow.profile\033[0m"
//     println " "
//     println "\033[2mCurrent User: $workflow.userName"
//     println "Nextflow-version: $nextflow.version"
//     println "Starting time: $nextflow.timestamp"
//     println "Workdir location:"
//     println "  $workflow.workDir\u001B[0m"
//     println " "

//     // error codes
//     if (params.profile) { exit 1, "--profile is WRONG use -profile" }
//     if ( !workflow.revision ) { 
//     println "\033[0;33mWARNING: It is recommended to use a stable release version via -r." 
//     println "Use 'nextflow info valegale/ONT_methylation' to check for available release versions.\033[0m\n"
//     }
//     // help
//     if (params.help) { exit 0, helpMSG() }

//     // --help
//     def helpMSG() {
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

// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/