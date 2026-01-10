#!/usr/bin/env python3
"""BAM → VCF germline calling with GATK HaplotypeCaller."""

from cosap_nx import BamReader, VariantCaller, Pipeline, PipelineRunner

BAM_PATH = "/home/ozan/grad_project/data/bam/NA12878_S1.chr20.10_10p1mb.bam"
REF_PATH = "/home/ozan/grad_project/data/refs/ucsc.hg19.chr20.unittest.fasta"
WORKDIR = "./workdir_hc"


def main():
    bam = BamReader(filename=BAM_PATH, name="NA12878_chr20")

    caller = VariantCaller(
        library="haplotypecaller",
        normal_sample=bam,
        params={"germline_sample_name": "NA12878_chr20"},
    )

    pipeline = Pipeline(ref_fasta=REF_PATH)
    pipeline.add(caller)

    config_path = pipeline.build(workdir=WORKDIR, cpus=4, memory="8 GB")

    runner = PipelineRunner()
    output_dir = runner.run_pipeline(config_path)

    print(f"\nDone. Output: {WORKDIR}/results/vcf/haplotypecaller/NA12878_chr20.vcf.gz")


if __name__ == "__main__":
    main()
