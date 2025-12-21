#!/usr/bin/env nextflow

/*
 * COSAP-NX v0.1 - Germline Variant Calling Pipeline
 * BAM -> DeepVariant -> VCF
 */

nextflow.enable.dsl = 2

// Validate required parameters
if (!params.bam) {
    error "Parameter 'bam' is required"
}
if (!params.ref_fasta) {
    error "Parameter 'ref_fasta' is required"
}

// Log pipeline info
log.info """
╔═══════════════════════════════════════════════════════════════╗
║                     COSAP-NX v0.1                             ║
║            Germline Variant Calling Pipeline                  ║
╚═══════════════════════════════════════════════════════════════╝

Parameters:
  - Sample ID    : ${params.sample_id}
  - BAM          : ${params.bam}
  - Reference    : ${params.ref_fasta}
  - Model type   : ${params.model_type}
  - Output GVCF  : ${params.gvcf}
  - Output dir   : ${params.outdir}
  - CPUs         : ${params.cpus}
  - Memory       : ${params.memory}
"""

/*
 * Process: DEEPVARIANT_CALL
 * Run DeepVariant for germline variant calling
 */
process DEEPVARIANT_CALL {
    tag "${sample_id}"
    // Using 1.6.0 - version 1.9.0 has a bug in postprocess_variants
    container 'google/deepvariant:1.6.0'
    cpus params.cpus
    memory params.memory

    publishDir "${params.outdir}/vcf/deepvariant", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bam_index)
    path ref_fasta
    path ref_fasta_fai
    path ref_fasta_dict

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: gvcf, optional: true

    script:
    def gvcf_arg = params.gvcf ? "--output_gvcf=${sample_id}.g.vcf.gz" : ""
    """
    # Use absolute path for intermediate results (like COSAP does)
    TMPDIR=\$(pwd)/tmpdir
    rm -rf \$TMPDIR
    mkdir -p \$TMPDIR

    run_deepvariant \\
        --model_type=${params.model_type} \\
        --ref=${ref_fasta} \\
        --reads=${bam} \\
        --output_vcf=${sample_id}.vcf.gz \\
        ${gvcf_arg} \\
        --num_shards=${task.cpus} \\
        --intermediate_results_dir=\$TMPDIR
    """
}

/*
 * Main workflow
 */
workflow {
    // Create input channels
    bam_ch = Channel.fromPath(params.bam, checkIfExists: true)

    // Find BAM index - try .bam.bai first, then .bai
    def bam_bai = file("${params.bam}.bai")
    def alt_bai = file("${params.bam}".replace('.bam', '.bai'))
    def bam_index_file = bam_bai.exists() ? bam_bai : alt_bai
    if (!bam_index_file.exists()) {
        error "BAM index not found. Tried: ${bam_bai} and ${alt_bai}"
    }
    bam_index_ch = Channel.fromPath(bam_index_file)

    ref_ch = Channel.fromPath(params.ref_fasta, checkIfExists: true)
    ref_fai_ch = Channel.fromPath("${params.ref_fasta}.fai", checkIfExists: true)

    // Try to find .dict file (optional for DeepVariant)
    def ref_dict_file = file("${params.ref_fasta}.dict")
    ref_dict_ch = ref_dict_file.exists() ? Channel.fromPath(ref_dict_file) : Channel.of(file('NO_DICT'))

    // Combine BAM with index
    bam_tuple = Channel.of(params.sample_id)
        .combine(bam_ch)
        .combine(bam_index_ch)

    // Run DeepVariant
    DEEPVARIANT_CALL(
        bam_tuple,
        ref_ch,
        ref_fai_ch,
        ref_dict_ch
    )

    // Log completion
    DEEPVARIANT_CALL.out.vcf.view { sample_id, vcf, tbi ->
        "Completed: ${sample_id} -> ${vcf}"
    }
}

workflow.onComplete {
    log.info """
Pipeline completed!
  - Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
  - Duration  : ${workflow.duration}
  - Output    : ${params.outdir}/vcf/deepvariant/
"""
}
