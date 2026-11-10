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

    ch_versions = channel.empty()

    // convert methylbam to fastq
    BAM2FASTQ(ch_samplesheet)
    fastq_ch = BAM2FASTQ.out.files
    ch_versions = ch_versions.mix(BAM2FASTQ.out.versions.first())
    
    // gzip reads
    // ZIPFASTQ(fastq_ch)

    // map fastq to reference assembly
    mapped_bams_ch = MINIMAP2(fastq_ch)

    // tabulate base modification calls across positions
    MODKIT_PILEUP(mapped_bams_ch)
    pileup_ch = MODKIT_PILEUP.out.files
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions.first())
    
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
    versions = ch_versions

}

workflow META_WORKFLOW {
    take:
    ch_samplesheet

    main:

    ch_versions = channel.empty()

    // convert methylbam to fastq
    BAM2FASTQ(ch_samplesheet)
    fastq_ch = BAM2FASTQ.out.files
    ch_versions = ch_versions.mix(BAM2FASTQ.out.versions.first())
    
    // gzip reads
    // ZIPFASTQ(fastq_ch)

    // map fastq to reference assembly
    mapped_bams_ch = MINIMAP2(fastq_ch)

    // Split bams by bin
    bins_ch = mapped_bams_ch.flatMap { meta, bam, index, assembly, bins ->

        def bin_files = bins.toFile().listFiles()
            ?.findAll { bin -> bin.name ==~ /.*\.(fasta|fa)$/ }
            ?.sort { bin -> bin.name }

        if (!bin_files) {
            error "No bin FASTA files found for sample ${meta.id}: ${bins}"
        }

        bin_files.collect { bin ->

        def bin_meta = meta + [
            bin_id: bin.baseName
        ]

        tuple(bin_meta, bam, index, assembly, file(bin))
        }
    }

    mapped_filtered_bams_ch = SPLIT_BAM_BY_BIN(bins_ch)

    // tabulate base modification calls across positions
    MODKIT_PILEUP(mapped_filtered_bams_ch)
    pileup_ch = MODKIT_PILEUP.out.files
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions.first())
    
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
    versions = ch_versions

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/