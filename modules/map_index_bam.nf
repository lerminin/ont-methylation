process BAM2FASTQ {
    // converts the BAM file from Dorado to FASTQ format 
    label 'minimap2'

    input:
    tuple val(meta), path(bam), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.fastq"), path(assembly)

    script:
    """
    samtools fastq ${bam} -T MM,ML > ${meta.id}.fastq
    """
    stub:
    """
    touch ${meta.id}.fastq    
    """
}

process ZIPFASTQ {
    // zip fastq file
    label 'minimap2'
    publishDir path: { "${params.outdir}/${meta.id}" }, mode:'copy'
        
    input:
    tuple val(meta), path(fastq_file), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), path(assembly)

    script:
    """
    gzip -c ${fastq_file} > ${meta.id}.fastq.gz
    """
    stub:
    """
    touch ${meta.id}.fastq.gz
    """
}

process MINIMAP2 {
    label 'minimap2'
    // align with minimap2 
    publishDir path: { "${params.outdir}/${meta.id}" }, mode:'copy'

    input:
    tuple val(meta), path(fastq_file), path(assembly)

    output:
    tuple val(meta), path("${meta.id}_methylation_mapped.bam"), path("${meta.id}_methylation_mapped.bam.bai"), path(assembly)

    script:
    """
    minimap2 -t ${task.cpus} --secondary=no -ax map-ont -y ${assembly} ${fastq_file} | \
    samtools view -b | \
    samtools sort -@ ${task.cpus} -o ${meta.id}_methylation_mapped.bam
    samtools index ${meta.id}_methylation_mapped.bam
    """
    stub:
    """
    touch ${meta.id}_methylation_mapped.bam
    touch ${meta.id}_methylation_mapped.bam.bai
    """
}


process SPLIT_BAM_BY_BIN {
    label 'minimap2'
    publishDir path: {"${params.outdir}/${meta.id}/bins/${bin_fa.baseName}"}, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bam_index), path(assembly), path(bin_fa)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), path(bin_fa)

    script:
    """
    # extract contigs from bin fasta
    awk '/^>/{gsub(/^>/,""); split(\$1,a," "); print a[1]}' ${bin_fa} > ${bin_fa.baseName}.contigs.txt

    # filter bam by contigs
    samtools view -b ${bam} \$(cat ${bin_fa.baseName}.contigs.txt) -o ${meta.id}.bam
    samtools index ${meta.id}.bam
    """
    stub:
    """
    touch ${meta.id}.bam
    touch ${meta.id}.bam.bai
    """
}
