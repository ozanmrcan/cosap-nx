#!/usr/bin/env python3
"""Alignment-only mode: FASTQ → BAM without variant calling."""

from cosap_nx import FastqReader, Mapper, Pipeline, PipelineRunner

FASTQ_R1 = "/home/ozan/grad_project/data/fastq/SRR1518011_1.fastq.gz"
FASTQ_R2 = "/home/ozan/grad_project/data/fastq/SRR1518011_2.fastq.gz"
REF_PATH = "/home/ozan/grad_project/data/refs/Homo_sapiens_assembly38.fasta"
WORKDIR = "/home/ozan/grad_project/cosap-nx/workdir_alignment"


def main():
    fastq_r1 = FastqReader(FASTQ_R1, read=1, name="NA12878_exome")
    fastq_r2 = FastqReader(FASTQ_R2, read=2, name="NA12878_exome")

    mapper = Mapper(
        library="bwa",
        input_step=[fastq_r1, fastq_r2],
        params={
            "read_groups": {
                "SM": "NA12878",
                "ID": "SRR1518011",
                "PL": "ILLUMINA",
                "LB": "exome",
            }
        },
    )

    pipeline = Pipeline(ref_fasta=REF_PATH)
    pipeline.add(mapper)

    config_path = pipeline.build(workdir=WORKDIR, cpus=4, memory="12 GB")

    runner = PipelineRunner()
    output_dir = runner.run_pipeline(config_path)

    print(f"\nDone. Output: {WORKDIR}/results/bam/NA12878_exome.bam")


if __name__ == "__main__":
    main()
