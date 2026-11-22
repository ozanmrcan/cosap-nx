/*
 * BWA-MEM alignment, sorted and indexed with samtools
 */
process BWA {
    tag "${meta.id}"
    // Mulled container with bwa 0.7.17 + samtools 1.15.1
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:66ed1b38d280722529bb8a0167b0cf02f8a0b488-0'
    cpus   { meta.cpus ?: params.defaults.cpus }
    memory { meta.memory ?: params.defaults.memory }
    publishDir "${params.outdir}/bam/bwa", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path fasta
    path ref_index

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai")

    script:
    // Read group from params.read_groups, like COSAP's Mapper
    def rg = [ID: meta.sample, SM: meta.sample, PL: 'ILLUMINA', LB: meta.sample] + (meta.params.read_groups ?: [:])
    def rg_line = "@RG\\t" + rg.collect { k, v -> "${k}:${v}" }.join("\\t")
    """
    bwa mem -t ${task.cpus} -R '${rg_line}' ${fasta} ${reads[0]} ${reads[1]} | \\
        samtools sort -@ ${task.cpus} -o ${meta.id}.bam -
    samtools index ${meta.id}.bam
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bam.bai
    """
}
