process MODKIT_PILEUP {
    label 'modkit'
    label 'process_high'
    publishDir {
        meta.bin_id
            ? "${params.outdir}/${meta.id}/${meta.bin_id}"
            : "${params.outdir}/${meta.id}"
    }, mode: 'copy', pattern: "*.bed"
    // execute the modkit pileup command
    
    input:
    tuple val(meta), path(mapped_bam), path(mapped_bam_bai), path(assembly)

    output:
    tuple val(meta), path("modkit_pileup_output.bed"), path(assembly), emit: files
    path "versions.yml", emit: versions

    script:
    filter_threshold = params.automatic_threshold_modkit ? '' : "--filter-threshold ${params.filter_threshold_modkit}"
    """
    modkit pileup -t ${task.cpus} ${mapped_bam} modkit_pileup_output.bed ${filter_threshold}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version | head -n 1)
    END_VERSIONS
    """
    stub:
    """
    touch modkit_pileup_output.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: "stub"
    END_VERSIONS   
    """
}

process MODKIT_PILEUP_BEDGRAPHS {
    label 'modkit'
    label 'process_high'
    publishDir {
        meta.bin_id
            ? "${params.outdir}/${meta.id}/${meta.bin_id}"
            : "${params.outdir}/${meta.id}"
    }, mode: 'copy'
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
    label 'process_low'
    publishDir {
        meta.bin_id
            ? "${params.outdir}/${meta.id}/${meta.bin_id}"
            : "${params.outdir}/${meta.id}"
    }, mode: 'copy'
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
    label 'modkit'
    label 'process_high'
    publishDir {
        meta.bin_id
            ? "${params.outdir}/${meta.id}/${meta.bin_id}"
            : "${params.outdir}/${meta.id}"
    }, mode: 'copy'
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

