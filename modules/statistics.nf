process COMPUTE_STATISTICS {
    label 'biopython'
    publishDir path: { "${params.outdir}/${meta.id}" }, mode:'copy'
    // run an inhouse script that computes which bases are methylated: modified bases / total bases 
    
    input:
    tuple val(meta), path(bed_file), path(assembly)

    output:
    tuple val(meta), path("methylation_statistics")

    script:
    """
    compute_statistics.py ${bed_file} ${assembly} methylation_statistics
    """
    stub:
    """
    mkdir -p methylation_statistics    
    """
}