#!/usr/bin/env python3
"""
COSAP-NX Example: Germline Variant Calling from FASTQ

This example demonstrates the complete FASTQ-to-VCF workflow:
1. Starting from paired-end FASTQ files
2. BWA-MEM alignment to reference genome
3. DeepVariant germline calling
4. Outputs: aligned BAM + variants VCF

Prerequisites:
- BWA index for reference genome (run: bwa index reference.fasta)
- Test FASTQ files (can be generated from BAM using samtools fastq)
"""

import os
from cosap_nx import FastqReader, VariantCaller, Pipeline, PipelineRunner


def main():
    # ====================================================================
    # Configuration - Update these paths for your data
    # ====================================================================

    # Input FASTQ files (paired-end)
    fastq_r1 = "/home/ozan/grad_project/data/fastq/NA12878_chr20_R1.fastq.gz"
    fastq_r2 = "/home/ozan/grad_project/data/fastq/NA12878_chr20_R2.fastq.gz"

    # Reference genome (must be BWA-indexed)
    ref_fasta = "/home/ozan/grad_project/data/refs/ucsc.hg19.chr20.unittest.fasta"


    # Working directory for outputs
    workdir = "/home/ozan/grad_project/cosap-nx/workdir_fastq"

    # Pipeline resources
    cpus = 8
    memory = "12 GB"  # Adjusted for systems with <16GB RAM

    # ====================================================================
    # Build and Run Pipeline
    # ====================================================================

    print("=" * 70)
    print("COSAP-NX v0.1.2: FASTQ → BAM → VCF Pipeline")
    print("=" * 70)

    # Create FASTQ readers for paired-end reads
    print(f"\nInput files:")
    print(f"  Read 1: {fastq_r1}")
    print(f"  Read 2: {fastq_r2}")

    fastq_1 = FastqReader(fastq_r1, read=1, name="NA12878_chr20")
    fastq_2 = FastqReader(fastq_r2, read=2, name="NA12878_chr20")

    # Configure variant caller with FASTQ input
    print(f"\nVariant Caller:")
    print(f"  Library: DeepVariant")
    print(f"  Model: WGS")
    print(f"  Resources: {cpus} CPUs, {memory} memory")

    caller = VariantCaller(
        library="deepvariant",
        normal_sample=[fastq_1, fastq_2],  # List of FastqReader instances
        params={"model_type": "WGS"},
        gvcf=False,
    )

    # Build pipeline
    print(f"\nBuilding pipeline...")

    pipeline = Pipeline(ref_fasta=ref_fasta)
    pipeline.add(caller)

    config_path = pipeline.build(workdir=workdir, cpus=cpus, memory=memory)

    print(f"✓ Config generated: {config_path}")

    # Run pipeline
    print(f"\n" + "=" * 70)
    print("Running Nextflow pipeline...")
    print("=" * 70)
    print("\nNote: This will take several minutes for alignment + variant calling.")
    print("Outputs:")
    print(f"  - Aligned BAM: {workdir}/results/bam/NA12878_chr20.bam")
    print(f"  - Variants VCF: {workdir}/results/vcf/deepvariant/NA12878_chr20.vcf.gz")
    print()

    runner = PipelineRunner()
    output_dir = runner.run_pipeline(config_path)

    print("\n" + "=" * 70)
    print("Pipeline Completed Successfully!")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir}")
    print(f"\nGenerated files:")
    print(f"  1. Aligned BAM: {output_dir}/bam/NA12878_chr20.bam")
    print(f"  2. BAM index: {output_dir}/bam/NA12878_chr20.bam.bai")
    print(f"  3. Variants VCF: {output_dir}/vcf/deepvariant/NA12878_chr20.vcf.gz")
    print(f"  4. VCF index: {output_dir}/vcf/deepvariant/NA12878_chr20.vcf.gz.tbi")
    print()


if __name__ == "__main__":
    main()
