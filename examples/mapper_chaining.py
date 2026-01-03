#!/usr/bin/env python3
"""
COSAP-NX Example: Mapper → VariantCaller Chaining

This example demonstrates the modular architecture by chaining
Mapper and VariantCaller components.

This shows how COSAP-NX components can be composed together,
similar to COSAP's design.
"""

from cosap_nx import FastqReader, Mapper, VariantCaller, Pipeline, PipelineRunner

# Input files
FASTQ_R1 = "/home/ozan/grad_project/data/fastq/SRR1518011_1.fastq.gz"
FASTQ_R2 = "/home/ozan/grad_project/data/fastq/SRR1518011_2.fastq.gz"
REF_PATH = "/home/ozan/grad_project/data/refs/Homo_sapiens_assembly38.fasta"
WORKDIR = "/home/ozan/grad_project/cosap-nx/workdir_chained"


def main():
    """Demonstrate Mapper → VariantCaller chaining."""

    # Step 1: Define FASTQ inputs
    fastq_r1 = FastqReader(FASTQ_R1, read=1, name="NA12878_exome")
    fastq_r2 = FastqReader(FASTQ_R2, read=2, name="NA12878_exome")

    # Step 2: Configure alignment
    mapper = Mapper(
        library="bwa",
        input_step=[fastq_r1, fastq_r2],
        params={
            "read_groups": {
                "SM": "NA12878",
                "ID": "SRR1518011",
                "PL": "ILLUMINA",
                "LB": "exome"
            }
        }
    )

    # Step 3: Configure variant calling using Mapper output
    # This demonstrates component chaining - Mapper feeds into VariantCaller
    variant_caller = VariantCaller(
        library="deepvariant",
        normal_sample=mapper,  # <-- Mapper as input (modular design!)
        params={
            "germline_sample_name": "NA12878_exome",
            "model_type": "WES"  # Whole Exome Sequencing
        },
        gvcf=True
    )

    # Step 4: Build pipeline
    # Note: We only add the VariantCaller - it internally uses Mapper
    # This is the COSAP-style chaining approach
    pipeline = Pipeline(ref_fasta=REF_PATH)
    pipeline.add(variant_caller)

    config_path = pipeline.build(
        workdir=WORKDIR,
        cpus=4,
        memory="14 GB"
    )

    # Step 5: Run the pipeline
    print("\n" + "="*60)
    print("Running Chained Pipeline: Mapper → VariantCaller")
    print("="*60)
    print("Architecture: Modular component composition")
    print(f"  1. Mapper (BWA)        : FASTQ → BAM")
    print(f"  2. VariantCaller (DV)  : BAM → VCF")
    print("="*60 + "\n")

    runner = PipelineRunner()
    output_dir = runner.run_pipeline(config_path)

    print("\n" + "="*60)
    print("Pipeline Complete!")
    print("="*60)
    print(f"Output directory: {output_dir}")
    print(f"VCF: {output_dir}/vcf/deepvariant/")
    print("="*60)


if __name__ == "__main__":
    main()
