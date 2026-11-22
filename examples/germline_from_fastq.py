#!/usr/bin/env python3
"""FASTQ → VCF germline calling pipeline."""

from cosap_nx import FastqReader, Mapper, VariantCaller, Pipeline, PipelineRunner


def main():
    # Inputs
    fastq_r1 = "/home/ozan/grad_project/data/fastq/NA12878_chr20_R1.fastq.gz"
    fastq_r2 = "/home/ozan/grad_project/data/fastq/NA12878_chr20_R2.fastq.gz"
    ref_fasta = "/home/ozan/grad_project/data/refs/ucsc.hg19.chr20.unittest.fasta"
    workdir = "/home/ozan/grad_project/cosap-nx/workdir_fastq"

    # Setup pipeline
    fastq_1 = FastqReader(fastq_r1, read=1, name="NA12878_chr20")
    fastq_2 = FastqReader(fastq_r2, read=2, name="NA12878_chr20")

    mapper = Mapper(
        library="bwa",
        input_step=[fastq_1, fastq_2],
        params={"read_groups": {"SM": "NA12878_chr20"}},
    )

    caller = VariantCaller(
        library="deepvariant",
        germline=mapper,
        params={"model_type": "WGS"},
    )

    pipeline = Pipeline(ref_fasta=ref_fasta)
    pipeline.add(caller)

    config_path = pipeline.build(workdir=workdir, cpus=8, memory="12 GB")

    # Run
    runner = PipelineRunner()
    output_dir = runner.run_pipeline(config_path)

    print(f"\nDone. Outputs:")
    print(f"  BAM: {workdir}/results/bam/bwa/NA12878_chr20_bwa.bam")
    print(f"  VCF: {workdir}/results/vcf/deepvariant/NA12878_chr20_bwa_deepvariant.vcf.gz")


if __name__ == "__main__":
    main()
