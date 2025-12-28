#!/usr/bin/env python3
"""
COSAP-NX Example: Homework Exome Data (NA12878_exome.bam)

This example tests COSAP-NX with the same data used in BLG348E homework:
  BAM (exome) -> DeepVariant (WES model) -> VCF + gVCF

Reference: Homo_sapiens_assembly38.fasta
Sample: NA12878_exome.bam

Usage:
    python homework_exome.py

Before running:
    1. Ensure files were copied to grad_project/data/
    2. Ensure Docker is running
    3. Ensure Nextflow is installed
"""

from cosap_nx import BamReader, VariantCaller, Pipeline, PipelineRunner

# ============================================================================
# Configuration - Homework Files
# ============================================================================

BAM_PATH = "/home/ozan/grad_project/data/bam/NA12878_exome.bam"
REF_PATH = "/home/ozan/grad_project/data/refs/Homo_sapiens_assembly38.fasta"
WORKDIR = "/home/ozan/grad_project/cosap-nx/workdir_homework"

# ============================================================================
# Pipeline Definition
# ============================================================================

def main():
    print("=" * 70)
    print("COSAP-NX: Homework Exome Test (Matching BLG348E COSAP Pipeline)")
    print("=" * 70)
    print()

    # 1. Define input BAM
    bam = BamReader(
        filename=BAM_PATH,
        name="NA12878_exome"
    )

    print(f"Input BAM: {BAM_PATH}")
    print(f"Reference: {REF_PATH}")
    print()

    # 2. Configure variant caller (matching homework COSAP settings)
    dv = VariantCaller(
        library="deepvariant",
        normal_sample=bam,
        params={
            "germline_sample_name": "NA12878_exome",
            "model_type": "WES",  # Whole Exome Sequencing (NOT WGS!)
        },
        gvcf=True  # Generate both VCF and gVCF (matching homework)
    )

    # 3. Build pipeline
    pipeline = Pipeline(ref_fasta=REF_PATH)
    pipeline.add(dv)

    # 4. Generate config
    config_path = pipeline.build(
        workdir=WORKDIR,
        cpus=8,
        memory="14 GB"
    )

    print(f"Pipeline config: {config_path}")
    print()

    # 5. Run the pipeline
    print("=" * 70)
    print("Running Nextflow pipeline...")
    print("=" * 70)
    print()
    print("Expected outputs:")
    print(f"  VCF:  {WORKDIR}/results/vcf/deepvariant/NA12878_exome.vcf.gz")
    print(f"  gVCF: {WORKDIR}/results/vcf/deepvariant/NA12878_exome.g.vcf.gz")
    print()

    runner = PipelineRunner()
    output_dir = runner.run_pipeline(config_path)

    print()
    print("=" * 70)
    print("Pipeline Complete!")
    print("=" * 70)
    print(f"Output directory: {output_dir}")
    print()
    print("Compare with COSAP homework results using:")
    print("  bcftools stats <vcf_file>")
    print("  bcftools isec <cosap_vcf> <cosapnx_vcf> -p comparison/")
    print()


if __name__ == "__main__":
    main()
