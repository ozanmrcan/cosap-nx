#!/usr/bin/env python3
"""
COSAP-NX Example: Germline Variant Calling with DeepVariant

This example demonstrates the v0.1 pipeline:
  BAM -> DeepVariant -> VCF

Usage:
    python germline_deepvariant.py

Before running:
    1. Set BAM_PATH to your aligned BAM file
    2. Set REF_PATH to your reference genome (with .fai index)
    3. Ensure Docker is running
    4. Ensure Nextflow is installed
"""

from cosap_nx import BamReader, VariantCaller, Pipeline, PipelineRunner

# ============================================================================
# Configuration - EDIT THESE PATHS
# ============================================================================

BAM_PATH = "/home/ozan/grad_project/data/bam/NA12878_S1.chr20.10_10p1mb.bam"
REF_PATH = "/home/ozan/grad_project/data/refs/ucsc.hg19.chr20.unittest.fasta"
WORKDIR = "./workdir"

# ============================================================================
# Pipeline Definition
# ============================================================================

def main():
    # 1. Define input BAM
    bam = BamReader(
        filename=BAM_PATH,
        name="NA12878_chr20"  # Optional: auto-derived from filename if omitted
    )

    # 2. Configure variant caller
    dv = VariantCaller(
        library="deepvariant",
        normal_sample=bam,
        params={
            "germline_sample_name": "NA12878_chr20",
            "model_type": "WGS",  # WGS for this test data
        },
        gvcf=False  # Set True to also output gVCF
    )

    # 3. Build pipeline
    pipeline = Pipeline()
    pipeline.set_reference(REF_PATH)
    pipeline.add(dv)

    # 4. Generate config
    config_path = pipeline.build(
        workdir=WORKDIR,
        cpus=4,
        memory="12 GB"
    )

    print(f"Config generated: {config_path}")
    print()

    # 5. Run the pipeline
    runner = PipelineRunner()
    output_dir = runner.run_pipeline(
        pipeline_config=config_path,
        profile="docker"
    )

    print()
    print(f"Pipeline complete! VCF output: {output_dir}")


if __name__ == "__main__":
    main()
