process compute_statistics {
    label 'biopython'
    // run an inhouse script that computes which bases are methylated: modified bases / total bases 
    
    input:
    tuple val(assembly_name), val(sample), path(assembly), path(bed_file)

    output:
    tuple val(assembly_name), path("methylation_statistics")

    script:
    """
    compute_statistics.py ${bed_file} ${assembly} methylation_statistics
    """
    stub:
    """
    mkdir -p methylation_statistics    
    """
}