process COMPUTE_STATISTICS {
    label 'biopython'
    publishDir {
        meta.bin_id
            ? "${params.outdir}/${meta.id}/${meta.bin_id}"
            : "${params.outdir}/${meta.id}"
    }, mode: 'copy'
    // run an inhouse script that computes which bases are methylated: modified bases / total bases 
    
    input:
    tuple val(meta), path(bed_file), path(assembly), path(bins)

    output:
    tuple val(meta), path("${prefix}_methylation_statistics")

    script:
    prefix = meta.bin_id ? "${meta.id}.${meta.bin_id}" : meta.id
    """
    compute_statistics.py ${bed_file} ${assembly} ${prefix}_methylation_statistics
    """
    stub:
    prefix = meta.bin_id ? "${meta.id}.${meta.bin_id}" : meta.id
    """
    mkdir -p ${prefix}_methylation_statistics    
    """
}