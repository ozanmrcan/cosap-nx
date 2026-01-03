#!/usr/bin/env nextflow

/*
 * COSAP-NX v0.1.2 - Germline Variant Calling Pipeline
 * Modes:
 *   1. FASTQ -> BWA-MEM -> sorted BAM (alignment only)
 *   2. FASTQ -> BWA-MEM -> sorted BAM -> DeepVariant -> VCF
 *   3. BAM -> DeepVariant -> VCF (backward compatible with v0.1.0)
 */

nextflow.enable.dsl = 2

// Validate required parameters
def alignment_only = params.mode == "alignment_only"

if (alignment_only) {
    // Alignment-only mode: requires FASTQ
    if (!params.fastq_r1 || !params.fastq_r2) {
        error "Alignment-only mode requires both --fastq_r1 and --fastq_r2"
    }
} else {
    // Variant calling mode: requires BAM or FASTQ
    if (!params.bam && (!params.fastq_r1 || !params.fastq_r2)) {
        error "Must provide either --bam or both --fastq_r1 and --fastq_r2"
    }
}

if (!params.ref_fasta) {
    error "Parameter 'ref_fasta' is required"
}

// Log pipeline info
def input_type = params.bam ? "BAM" : "FASTQ"
def input_files = params.bam ? params.bam : "${params.fastq_r1}, ${params.fastq_r2}"

log.info """
╔═══════════════════════════════════════════════════════════════╗
║                    COSAP-NX v0.1.2                            ║
║            Germline Variant Calling Pipeline                  ║
╚═══════════════════════════════════════════════════════════════╝

Parameters:
  - Sample ID    : ${params.sample_id}
  - Input type   : ${input_type}
  - Input files  : ${input_files}
  - Reference    : ${params.ref_fasta}
  - Model type   : ${params.model_type}
  - Output GVCF  : ${params.gvcf}
  - Output dir   : ${params.outdir}
  - CPUs         : ${params.cpus}
  - Memory       : ${params.memory}
"""

/*
 * Process: BWA_ALIGN
 * Align FASTQ reads to reference genome using BWA-MEM
 */
process BWA_ALIGN {
    tag "${sample_id}"
    // Mulled container with bwa 0.7.17 + samtools 1.15.1
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:66ed1b38d280722529bb8a0167b0cf02f8a0b488-0'
    cpus params.cpus
    memory params.memory

    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path ref_fasta
    path ref_fasta_fai
    path bwa_index_files  // .amb, .ann, .bwt, .pac, .sa (as list)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    script:
    // Use custom read group if provided (for Mapper), otherwise use default
    def read_group = params.read_group ?: "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}"
    """
    # Align with BWA-MEM, add read groups, pipe to samtools sort
    bwa mem -t ${task.cpus} -R '${read_group}' \\
        ${ref_fasta} ${read1} ${read2} | \\
    samtools sort -@ ${task.cpus} -o ${sample_id}.bam -

    # Index the sorted BAM
    samtools index ${sample_id}.bam
    """
}

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
    // Common reference channels
    ref_ch = Channel.fromPath(params.ref_fasta, checkIfExists: true)
    ref_fai_ch = Channel.fromPath("${params.ref_fasta}.fai", checkIfExists: true)

    // Try to find .dict file (optional for DeepVariant)
    def ref_dict_file = file("${params.ref_fasta}.dict")
    ref_dict_ch = ref_dict_file.exists() ? Channel.fromPath(ref_dict_file) : Channel.of(file('NO_DICT'))

    // Branch based on input type and mode (using alignment_only defined at top)
    if (alignment_only) {
        // Mode 1: Alignment only (FASTQ -> BAM)
        log.info "Running Alignment-Only mode: FASTQ -> BAM"

        // Check for BWA index files (also checks .64 suffix for large genomes)
        bwa_index_files = Channel.fromPath("${params.ref_fasta}.{amb,ann,bwt,pac,sa,64.amb,64.ann,64.bwt,64.pac,64.sa}")
            .collect()
            .ifEmpty { error "BWA index not found. Run: bwa index ${params.ref_fasta}" }

        // Create FASTQ input channels
        fastq_r1_ch = Channel.fromPath(params.fastq_r1, checkIfExists: true)
        fastq_r2_ch = Channel.fromPath(params.fastq_r2, checkIfExists: true)

        // Combine into tuple (sample_id, read1, read2)
        fastq_tuple = Channel.of(params.sample_id)
            .combine(fastq_r1_ch)
            .combine(fastq_r2_ch)

        // Run alignment ONLY (no variant calling)
        BWA_ALIGN(fastq_tuple, ref_ch, ref_fai_ch, bwa_index_files)

        // Log completion
        BWA_ALIGN.out.bam.view { sample_id, bam, bai ->
            "Completed: ${sample_id} -> ${bam}"
        }

    } else if (params.fastq_r1 && params.fastq_r2) {
        // Mode 2: Full pipeline (FASTQ -> BAM -> VCF)
        log.info "Running FASTQ mode: Alignment + Variant Calling"

        // Check for BWA index files (also checks .64 suffix for large genomes)
        bwa_index_files = Channel.fromPath("${params.ref_fasta}.{amb,ann,bwt,pac,sa,64.amb,64.ann,64.bwt,64.pac,64.sa}")
            .collect()
            .ifEmpty { error "BWA index not found. Run: bwa index ${params.ref_fasta}" }

        // Create FASTQ input channels
        fastq_r1_ch = Channel.fromPath(params.fastq_r1, checkIfExists: true)
        fastq_r2_ch = Channel.fromPath(params.fastq_r2, checkIfExists: true)

        // Combine into tuple (sample_id, read1, read2)
        fastq_tuple = Channel.of(params.sample_id)
            .combine(fastq_r1_ch)
            .combine(fastq_r2_ch)

        // Run alignment
        BWA_ALIGN(fastq_tuple, ref_ch, ref_fai_ch, bwa_index_files)

        // Feed aligned BAM to variant caller
        DEEPVARIANT_CALL(BWA_ALIGN.out.bam, ref_ch, ref_fai_ch, ref_dict_ch)

        // Log completion
        DEEPVARIANT_CALL.out.vcf.view { sample_id, vcf, tbi ->
            "Completed: ${sample_id} -> ${vcf}"
        }

    } else if (params.bam) {
        // Mode 3: BAM mode (BAM -> VCF, backward compatible with v0.1.0)
        log.info "Running BAM mode: Variant Calling only"

        // Create BAM input channels
        bam_ch = Channel.fromPath(params.bam, checkIfExists: true)

        // Find BAM index - try .bam.bai first, then .bai
        def bam_bai = file("${params.bam}.bai")
        def alt_bai = file("${params.bam}".replace('.bam', '.bai'))
        def bam_index_file = bam_bai.exists() ? bam_bai : alt_bai
        if (!bam_index_file.exists()) {
            error "BAM index not found. Tried: ${bam_bai} and ${alt_bai}"
        }
        bam_index_ch = Channel.fromPath(bam_index_file)

        // Combine BAM with index
        bam_tuple = Channel.of(params.sample_id)
            .combine(bam_ch)
            .combine(bam_index_ch)

        // Run DeepVariant
        DEEPVARIANT_CALL(bam_tuple, ref_ch, ref_fai_ch, ref_dict_ch)

        // Log completion
        DEEPVARIANT_CALL.out.vcf.view { sample_id, vcf, tbi ->
            "Completed: ${sample_id} -> ${vcf}"
        }
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
