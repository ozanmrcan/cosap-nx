/*
 * GATK HaplotypeCaller germline calling.
 * With gvcf, only the gVCF is written (same as COSAP).
 */
process HAPLOTYPECALLER {
    tag "${meta.id}"
    container 'broadinstitute/gatk:4.5.0.0'
    cpus   { meta.cpus ?: params.defaults.cpus }
    memory { meta.memory ?: params.defaults.memory }
    publishDir "${params.outdir}/vcf/haplotypecaller", mode: 'copy'

    input:
    tuple val(meta), path(bam_files)
    path fasta
    path ref_index

    output:
    tuple val(meta), path("${meta.id}.{vcf,g.vcf}.gz"), path("${meta.id}.{vcf,g.vcf}.gz.tbi")

    script:
    def bam = bam_files.find { it.name.endsWith('.bam') }
    def out = meta.gvcf ? "${meta.id}.g.vcf.gz" : "${meta.id}.vcf.gz"
    def erc = meta.gvcf ? "-ERC GVCF" : ""
    """
    # GATK expects ref.dict next to ref.fasta, link it if only ref.fasta.dict exists
    DICT=\$(basename ${fasta} .fasta).dict
    [ -e \$DICT ] || ln -s ${fasta}.dict \$DICT

    gatk HaplotypeCaller \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${out} \\
        ${erc} \\
        --native-pair-hmm-threads ${task.cpus}
    """

    stub:
    """
    touch ${meta.id}.vcf.gz ${meta.id}.vcf.gz.tbi
    """
}
