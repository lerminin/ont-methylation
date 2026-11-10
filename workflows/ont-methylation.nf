/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM2FASTQ; ZIPFASTQ; MINIMAP2; SPLIT_BAM_BY_BIN } from '../modules/map_index_bam.nf'
include { MODKIT_PILEUP; MODKIT_PILEUP_BEDGRAPHS; MODKIT_FIND_MOTIFS; CUSTOM_BEDGRAPHS} from '../modules/modkit.nf'
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

workflow ISOLATE_WORKFLOW {
    take:
    ch_samplesheet

    main:

    // convert methylbam to fastq
    fastq_ch = BAM2FASTQ(ch_samplesheet)
    
    // gzip reads
    // ZIPFASTQ(fastq_ch)

    // map fastq to reference assembly
    mapped_bams_ch = MINIMAP2(fastq_ch)

    // tabulate base modification calls across positions
    pileup_ch = MODKIT_PILEUP(mapped_bams_ch)
    
    // generate bedgraphs
    pileup_bedgraphs_ch = MODKIT_PILEUP_BEDGRAPHS(mapped_bams_ch)
    
    // run script that computes which bases are methylated
    custom_bedgraphs_ch = CUSTOM_BEDGRAPHS(pileup_ch)

    // find methylated motifs
    motifs_ch = MODKIT_FIND_MOTIFS(pileup_ch)

    // count statistics and percent bases modified
    statistics_ch = COMPUTE_STATISTICS(pileup_ch)

    emit:
    pileup = pileup_ch
    bedgraphs = pileup_bedgraphs_ch
    custom_bedgraphs = custom_bedgraphs_ch
    statistics = statistics_ch
    motifs = motifs_ch

}

workflow META_WORKFLOW {
    take:
    ch_samplesheet

    main:

    // convert methylbam to fastq
    fastq_ch = BAM2FASTQ(ch_samplesheet)
    
    // gzip reads
    // ZIPFASTQ(fastq_ch)

    // map fastq to reference assembly
    mapped_bams_ch = MINIMAP2(fastq_ch)

    // Split bams by bin
    bins_ch = mapped_bams_ch.map { meta, _bam, _index, _assembly, bins ->
            tuple(meta.id, bins)
        }

    // Match mapped BAMs to the correct sample's bins
    mapped_keyed_ch = mapped_bams_ch.map { meta, bam, index, assembly ->
        tuple(meta.id, meta, bam, index, assembly)
    }

    bam_bin_pairs = mapped_keyed_ch.join(bins_ch)
    mapped_filtered_bams_ch = SPLIT_BAM_BY_BIN(bam_bin_pairs)

    // tabulate base modification calls across positions
    pileup_ch = MODKIT_PILEUP(mapped_filtered_bams_ch)
    
    // generate bedgraphs
    pileup_bedgraphs_ch = MODKIT_PILEUP_BEDGRAPHS(mapped_filtered_bams_ch)
    
    // run script that computes which bases are methylated
    custom_bedgraphs_ch = CUSTOM_BEDGRAPHS(pileup_ch)

    // find methylated motifs
    motifs_ch = MODKIT_FIND_MOTIFS(pileup_ch)

    // count statistics and percent bases modified
    statistics_ch = COMPUTE_STATISTICS(pileup_ch)

    emit:
    pileup = pileup_ch
    bedgraphs = pileup_bedgraphs_ch
    custom_bedgraphs = custom_bedgraphs_ch
    statistics = statistics_ch
    motifs = motifs_ch

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/