process bam2fastq {
    // converts the BAM file from Dorado to FASTQ format 
    label 'minimap2'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.fastq")

    script:
    """
    samtools fastq ${bam} -T MM,ML > ${sample}.fastq
    """
    stub:
    """
    touch ${sample}.fastq    
    """
}

process zipfastq {
    // zip fastq file
    label 'minimap2'
    publishDir  "${params.outdir}/${sample}", mode:'copy'
        
    input:
    tuple val(sample), path(fastq_file)

    output:
    path "${sample}.fastq.gz"

    script:
    """
    gzip -c ${fastq_file} > ${sample}.fastq.gz
    """
    stub:
    """
    touch ${sample}.fastq.gz
    """
}

process minimap2 {
    label 'minimap2'
    // align with minimap2 
    publishDir  "${params.outdir}/${sample}", mode:'copy'

    input:
    tuple val(sample), path(fastq_file), path(assembly)

    output:
    tuple val(sample), path("${sample}_methylation_mapped.bam"), path("${sample}_methylation_mapped.bam.bai"), path(assembly)

    script:
    """
    minimap2 -t ${task.cpus} --secondary=no -ax map-ont -y ${assembly} ${fastq_file} | \
    samtools view -b | \
    samtools sort -@ ${task.cpus} -o ${sample}_methylation_mapped.bam
    samtools index ${sample}_methylation_mapped.bam
    """
    stub:
    """
    touch ${sample}_methylation_mapped.bam
    touch ${sample}_methylation_mapped.bam.bai
    """
}


process split_bam_by_bin {
    label 'minimap2'
    publishDir "${params.outdir}/${sample}/bins/${bin_fa.baseName}", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bam_index), path(assembly), path(bin_fa)

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"), path(bin_fa)

    script:
    """
    # extract contigs from bin fasta
    awk '/^>/{gsub(/^>/,""); split(\$1,a," "); print a[1]}' ${bin_fa} > ${bin_fa.baseName}.contigs.txt

    # filter bam by contigs
    samtools view -b ${bam} \$(cat ${bin_fa.baseName}.contigs.txt) -o ${sample}.bam
    samtools index ${sample}.bam
    """
    stub:
    """
    touch ${sample}.bam
    touch ${sample}.bam.bai
    """
}
