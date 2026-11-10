process MODKIT_PILEUP {
    label 'modkit'
    label 'process_high'
    publishDir {
        meta.bin_id
            ? "${params.outdir}/${meta.id}/${meta.bin_id}"
            : "${params.outdir}/${meta.id}"
    }, mode: 'copy'
    // execute the modkit pileup command
    
    input:
    tuple val(meta), path(mapped_bam), path(mapped_bam_bai), path(assembly), path(bins)

    output:
    tuple val(meta), path("${prefix}_modkit_pileup_output.bed"), path(assembly), path(bins), emit: files
    path "versions.yml", emit: versions

    script:
    filter_threshold = params.automatic_threshold_modkit ? '' : "--filter-threshold ${params.filter_threshold_modkit}"
    prefix = meta.bin_id ? "${meta.id}.${meta.bin_id}" : meta.id
    """
    modkit pileup -t ${task.cpus} ${mapped_bam} ${prefix}_modkit_pileup_output.bed ${filter_threshold}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version | head -n 1)
    END_VERSIONS
    """
    stub:
    prefix = meta.bin_id ? "${meta.id}.${meta.bin_id}" : meta.id
    """
    touch ${prefix}_modkit_pileup_output.bed

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
    tuple val(meta), path(mapped_bam), path(mapped_bam_bai), path(assembly), path(bins)
    output:
    tuple val(meta), path("bedgraphs_${prefix}"), path(assembly)

    script:
    prefix = meta.bin_id ? "${meta.id}.${meta.bin_id}" : meta.id
    """
    modkit pileup -t ${task.cpus} ${mapped_bam} --bedgraph bedgraphs_${prefix} --filter-threshold ${params.filter_threshold_modkit} 
    """
    stub:
    prefix = meta.bin_id ? "${meta.id}.${meta.bin_id}" : meta.id
    """
    mkdir -p bedgraphs_${prefix}
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
    tuple val(meta), path(bed_file), path(assembly), path(bins)

    output:
    tuple val(meta), path("bedgraphs_customized_${prefix}"), path("modifications_tables_${prefix}")

    script:
    prefix = meta.bin_id ? "${meta.id}.${meta.bin_id}" : meta.id
    """
    mkdir -p bedgraphs_customized_${prefix}
    mkdir -p modifications_tables_${prefix}

    custom_bedgraphs.py ${bed_file} ${assembly} . --percent_cutoff ${params.percent_cutoff_modification_table}
    """
    stub:
    prefix = meta.bin_id ? "${meta.id}.${meta.bin_id}" : meta.id
    """
    mkdir -p bedgraphs_customized_${prefix}
    mkdir -p modifications_tables_${prefix}
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
    tuple val(meta), path(bed_file), path(assembly), path(bins)

    output:
    tuple val(meta), path("${prefix}_modkit_motifs.tsv")

    script:
    prefix = meta.bin_id ? "${meta.id}.${meta.bin_id}" : meta.id
    """
    modkit find-motifs -t ${task.cpus} --in-bedmethyl ${bed_file} --ref ${assembly} -o ${prefix}_modkit_motifs.tsv
    """
    stub:
    prefix = meta.bin_id ? "${meta.id}.${meta.bin_id}" : meta.id
    """
    touch ${prefix}_modkit_motifs.tsv
    """ 
}

