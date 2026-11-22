/*
 * BWA-MEM2 alignment, sorted and indexed with samtools
 */
process BWA2 {
    tag "${meta.id}"
    // Mulled container with bwa-mem2 2.2.1 + samtools 1.19.2
    container 'quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2d15960ccea84e249a150b7f5d4db3a42fc2d6c3-0'
    cpus   { meta.cpus ?: params.defaults.cpus }
    memory { meta.memory ?: params.defaults.memory }
    publishDir "${params.outdir}/bam/bwa2", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path fasta
    path ref_index

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai")

    script:
    // Same read group handling as the BWA module
    def rg = [ID: meta.sample, SM: meta.sample, PL: 'ILLUMINA', LB: meta.sample] + (meta.params.read_groups ?: [:])
    def rg_line = "@RG\\t" + rg.collect { k, v -> "${k}:${v}" }.join("\\t")
    """
    bwa-mem2 mem -t ${task.cpus} -R '${rg_line}' ${fasta} ${reads[0]} ${reads[1]} | \\
        samtools sort -@ ${task.cpus} -o ${meta.id}.bam -
    samtools index ${meta.id}.bam
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bam.bai
    """
}
