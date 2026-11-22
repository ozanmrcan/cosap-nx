/*
 * DeepVariant germline calling
 */
process DEEPVARIANT {
    tag "${meta.id}"
    // Using 1.6.0 - version 1.9.0 has a bug in postprocess_variants
    container 'google/deepvariant:1.6.0'
    cpus   { meta.cpus ?: params.defaults.cpus }
    memory { meta.memory ?: params.defaults.memory }
    publishDir "${params.outdir}/vcf/deepvariant", mode: 'copy'

    input:
    tuple val(meta), path(bam_files)
    path fasta
    path ref_index

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), path("${meta.id}.vcf.gz.tbi")
    tuple val(meta), path("${meta.id}.g.vcf.gz"), path("${meta.id}.g.vcf.gz.tbi"), optional: true

    script:
    def bam = bam_files.find { it.name.endsWith('.bam') }
    def model_type = meta.params.model_type ?: 'WGS'
    def gvcf_arg = meta.gvcf ? "--output_gvcf=${meta.id}.g.vcf.gz" : ""
    """
    # Use absolute path for intermediate results (like COSAP does)
    TMPDIR=\$(pwd)/tmpdir
    rm -rf \$TMPDIR
    mkdir -p \$TMPDIR

    run_deepvariant \\
        --model_type=${model_type} \\
        --ref=${fasta} \\
        --reads=${bam} \\
        --output_vcf=${meta.id}.vcf.gz \\
        ${gvcf_arg} \\
        --num_shards=${task.cpus} \\
        --intermediate_results_dir=\$TMPDIR
    """

    stub:
    """
    touch ${meta.id}.vcf.gz ${meta.id}.vcf.gz.tbi
    """
}
