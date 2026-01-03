#!/usr/bin/env python3
"""
COSAP-NX Example: Alignment-Only Mode

This example demonstrates using the Mapper class for FASTQ → BAM alignment
without variant calling.

Use case: When you only need aligned BAM files and want to run variant
calling separately or with a different tool.
"""

from cosap_nx import FastqReader, Mapper, Pipeline, PipelineRunner

# Input files (homework exome FASTQ data)
FASTQ_R1 = "/home/ozan/grad_project/data/fastq/SRR1518011_1.fastq.gz"
FASTQ_R2 = "/home/ozan/grad_project/data/fastq/SRR1518011_2.fastq.gz"
REF_PATH = "/home/ozan/grad_project/data/refs/Homo_sapiens_assembly38.fasta"
WORKDIR = "/home/ozan/grad_project/cosap-nx/workdir_alignment"


def main():
    """Run alignment-only pipeline: FASTQ → BAM."""

    # Step 1: Define FASTQ inputs
    fastq_r1 = FastqReader(FASTQ_R1, read=1, name="NA12878_exome")
    fastq_r2 = FastqReader(FASTQ_R2, read=2, name="NA12878_exome")

    # Step 2: Configure alignment with BWA
    mapper = Mapper(
        library="bwa",
        input_step=[fastq_r1, fastq_r2],
        params={
            "read_groups": {
                "SM": "NA12878",        # Sample name (required)
                "ID": "SRR1518011",     # Run ID
                "PL": "ILLUMINA",       # Platform
                "LB": "exome"           # Library
            }
        }
    )

    # Step 3: Build pipeline (alignment only, no variant calling)
    pipeline = Pipeline(ref_fasta=REF_PATH)
    pipeline.add(mapper)

    config_path = pipeline.build(
        workdir=WORKDIR,
        cpus=4,
        memory="12 GB"
    )

    # Step 4: Run the pipeline
    print("\n" + "="*60)
    print("Running Alignment-Only Pipeline")
    print("="*60)
    print(f"Input:  {FASTQ_R1}")
    print(f"        {FASTQ_R2}")
    print(f"Output: {WORKDIR}/results/bam/NA12878_exome.bam")
    print("="*60 + "\n")

    runner = PipelineRunner()
    output_dir = runner.run_pipeline(config_path)

    print("\n" + "="*60)
    print("Alignment Complete!")
    print("="*60)
    print(f"Output directory: {output_dir}")
    print(f"BAM file: {output_dir}/bam/NA12878_exome.bam")
    print(f"BAM index: {output_dir}/bam/NA12878_exome.bam.bai")
    print("="*60)


if __name__ == "__main__":
    main()
