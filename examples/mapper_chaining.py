#!/usr/bin/env python3
"""Mapper → VariantCaller chaining example."""

from cosap_nx import FastqReader, Mapper, VariantCaller, Pipeline, PipelineRunner

FASTQ_R1 = "/home/ozan/grad_project/data/fastq/SRR1518011_1.fastq.gz"
FASTQ_R2 = "/home/ozan/grad_project/data/fastq/SRR1518011_2.fastq.gz"
REF_PATH = "/home/ozan/grad_project/data/refs/Homo_sapiens_assembly38.fasta"
WORKDIR = "/home/ozan/grad_project/cosap-nx/workdir_chained"


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

    # Mapper as input to VariantCaller (component chaining)
    caller = VariantCaller(
        library="deepvariant",
        normal_sample=mapper,
        params={"germline_sample_name": "NA12878_exome", "model_type": "WES"},
        gvcf=True,
    )

    pipeline = Pipeline(ref_fasta=REF_PATH)
    pipeline.add(caller)

    config_path = pipeline.build(workdir=WORKDIR, cpus=4, memory="14 GB")

    runner = PipelineRunner()
    output_dir = runner.run_pipeline(config_path)

    print(f"\nDone. Outputs:")
    print(f"  BAM:  {WORKDIR}/results/bam/NA12878_exome.bam")
    print(f"  VCF:  {WORKDIR}/results/vcf/deepvariant/NA12878_exome.vcf.gz")
    print(f"  gVCF: {WORKDIR}/results/vcf/deepvariant/NA12878_exome.g.vcf.gz")


if __name__ == "__main__":
    main()
