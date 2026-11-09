process MODKIT_PILEUP {
    label 'modkit'
    // execute the modkit pileup command
    
    input:
    tuple val(meta), path(mapped_bam), path(mapped_bam_bai), path(assembly)

    output:
    tuple val(meta), path("modkit_pileup_output.bed"), path(assembly)

    script:
    filter_threshold = params.automatic_threshold_modkit ? '' : "--filter-threshold ${params.filter_threshold_modkit}"
    """
    modkit pileup -t ${task.cpus} ${mapped_bam} modkit_pileup_output.bed ${filter_threshold}
    """
    stub:
    """
    touch modkit_pileup_output.bed   
    """
}

process MODKIT_PILEUP_BEDGRAPHS {
    label 'modkit'
    // execute the modkit pileup command to obtain the bedgraphs
    
    input:
    tuple val(meta), path(mapped_bam), path(mapped_bam_bai), path(assembly)

    output:
    tuple val(meta), path("bedgraphs"), path(assembly)

    script:
    """
    modkit pileup -t ${task.cpus} ${mapped_bam} --bedgraph bedgraphs --filter-threshold ${params.filter_threshold_modkit} 
    """
    stub:
    """
    mkdir -p bedgraphs
    """
}

process CUSTOM_BEDGRAPHS {
    label 'biopython'
    // run an inhouse script that computes which bases are methylated: modified bases / total bases.
    // additionally, it saves the positions with a high methylation levels (>0.5) in tsv tables, one for each modification (6mA, 5mC and 4mC).

    input:
    tuple val(meta), path(bed_file), path(assembly)

    output:
    tuple val(meta), path("bedgraphs_customized"), path("modifications_tables")

    script:
    """
    mkdir -p bedgraphs_customized
    mkdir -p modifications_tables

    custom_bedgraphs.py ${bed_file} ${assembly} . --percent_cutoff ${params.percent_cutoff_modification_table}
    """
    stub:
    """
    mkdir -p bedgraphs_customized
    mkdir -p modifications_tables
    """
}

process MODKIT_FIND_MOTIFS {
    label 'modkit_low'
    // find motifs from the output of modkit pileup
    //errorStrategy 'ignore'

    input:
    tuple val(meta), path(bed_file), path(assembly)

    output:
    tuple val(meta), path("modkit_motifs.tsv")

    script:
    """
    modkit find-motifs -t ${task.cpus} --in-bedmethyl ${bed_file} --ref ${assembly} -o modkit_motifs.tsv
    """
    stub:
    """
    touch modkit_motifs.tsv
    """ 
}


process PUBLISH_RESULTS_META {
    label 'publish'
    publishDir path: {"${params.outdir}/${meta.id}/bins/${assembly_name}"}, mode: 'copy'

    input:
    tuple val(assembly_name), val(meta), path(assembly), 
          path(bed_file), path(pileup_bedgraphs), path(custom_bedgraphs), path(modifications_tables), path(statistics)

    output:
    tuple path(bed_file), path(pileup_bedgraphs),  path(custom_bedgraphs), path(modifications_tables), path(statistics)

    script:
    """
    """
}

process PUBLISH_RESULTS_MOTIFS_META {
    label 'publish'
    publishDir path: {"${params.outdir}/${meta.id}/bins/${assembly_name}"}, mode: 'copy'

    input:
    tuple val(assembly_name), val(meta), path(assembly), 
          path(bed_file), path(motifs)

    output:
    path(motifs)

    script:
    """
    """
}

process PUBLISH_RESULTS {
    label 'publish'
    publishDir path: {"${params.outdir}/${meta.id}"}, mode: 'copy'

    input:
    tuple val(assembly_name), val(meta), path(assembly), 
          path(bed_file), path(pileup_bedgraphs), path(custom_bedgraphs), path(modifications_tables), path(statistics)

    output:
    tuple path(bed_file), path(pileup_bedgraphs), path(custom_bedgraphs), path(modifications_tables), path(statistics)

    script:
    """
    """
}

process PUBLISH_RESULTS_MOTIFS {
    label 'publish'
    publishDir path: {"${params.outdir}/${meta.id}"}, mode: 'copy'

    input:
    tuple val(assembly_name), val(meta), path(assembly), 
          path(bed_file), path(motifs)

    output:
    path(motifs)

    script:
    """
    """
}