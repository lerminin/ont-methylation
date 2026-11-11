process BAM2FASTQ {
    // converts the BAM file from Dorado to FASTQ format 
    label 'minimap2'
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.fastq"), path(assembly), emit: files
    path "versions.yml", emit: versions

    script:
    """
    samtools fastq ${bam} -T MM,ML > ${meta.id}.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1)
        minimap2: \$(minimap2 --version | head -n 1)
    END_VERSIONS
    """
    stub:
    """
    touch ${meta.id}.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "stub"
        minimap2: "stub"
    END_VERSIONS    
    """
}

process ZIPFASTQ {
    // zip fastq file
    label 'minimap2'
    label 'process_low'
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
    gzip -c /dev/null > ${meta.id}.fastq.gz
    """
}

process MINIMAP2 {
    label 'minimap2'
    label 'process_medium'
    // align with minimap2 
    publishDir path: { "${params.outdir}/${meta.id}" }, mode:'copy'

    input:
    tuple val(meta), path(fastq_file), path(assembly)

    output:
    tuple val(meta), path("methylation_mapped.bam"), path("methylation_mapped.bam.bai"), path(assembly)

    script:

    """
    minimap2 -t ${task.cpus} --secondary=no -ax map-ont -y ${assembly} ${fastq_file} | \
    samtools view -b | \
    samtools sort -@ ${task.cpus} -o methylation_mapped.bam
    samtools index methylation_mapped.bam
    """
    stub:
    """
    touch methylation_mapped.bam
    touch methylation_mapped.bam.bai
    """
}


process SPLIT_BAM_BY_BIN {
    label 'minimap2'
    label 'process_low'
    publishDir path: {"${params.outdir}/${meta.id}/${meta.bin_id}"}, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bam_index), path(assembly), path(bin_fasta)

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), path(bin_fasta)

    script:
    prefix = "${meta.id}.${meta.bin_id}"
    """
    # extract contigs from bin fasta
    awk '/^>/{gsub(/^>/,""); split(\$1,a," "); print a[1]}' ${bin_fasta} > ${prefix}.contigs.txt

    # filter bam by contigs
    samtools view -b ${bam} \$(cat ${prefix}.contigs.txt) -o ${prefix}.bam
    samtools index ${prefix}.bam
    """
    stub:
    prefix = "${meta.id}.${meta.bin_id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    """
}
