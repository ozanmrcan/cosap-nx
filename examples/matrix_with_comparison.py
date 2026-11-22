#!/usr/bin/env python3
"""FASTQ → {BWA, BWA-MEM2} → {DeepVariant, HaplotypeCaller} in one run, compared at the end."""

from cosap_nx import FastqReader, Mapper, VariantCaller, Pipeline, PipelineRunner


def main():
    # Inputs
    fastq_r1 = "/home/ozan/grad_project/data/fastq/NA12878_chr20_R1.fastq.gz"
    fastq_r2 = "/home/ozan/grad_project/data/fastq/NA12878_chr20_R2.fastq.gz"
    ref_fasta = "/home/ozan/grad_project/data/refs/ucsc.hg19.chr20.unittest.fasta"
    truth_vcf = "/home/ozan/grad_project/cosap-nx/tests/data/truth_chr20_10-10.1mb.vcf.gz"
    workdir = "/home/ozan/grad_project/cosap-nx/workdir_matrix"

    reads = [
        FastqReader(fastq_r1, read=1, name="NA12878_chr20"),
        FastqReader(fastq_r2, read=2, name="NA12878_chr20"),
    ]

    # 2 mappers x 2 callers
    pipeline = Pipeline(ref_fasta=ref_fasta)
    for mapper in ["bwa", "bwa2"]:
        mapped = Mapper(library=mapper, input_step=reads, params={"read_groups": {"SM": "NA12878_chr20"}})
        for caller in ["deepvariant", "haplotypecaller"]:
            pipeline.add(VariantCaller(library=caller, germline=mapped))

    config_path = pipeline.build(workdir=workdir, cpus=8, memory="12 GB")

    # Run, then compare the four VCFs
    results = PipelineRunner().run_pipeline(config_path, compare=True, truth_vcf=truth_vcf)

    print(f"\nDone. Comparison: {results}/comparison/NA12878_chr20/")


if __name__ == "__main__":
    main()
