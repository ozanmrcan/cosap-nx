#!/usr/bin/env python3
"""BLG348E homework exome data: BAM → VCF with WES model."""

from cosap_nx import BamReader, VariantCaller, Pipeline, PipelineRunner

BAM_PATH = "/home/ozan/grad_project/data/bam/NA12878_exome.bam"
REF_PATH = "/home/ozan/grad_project/data/refs/Homo_sapiens_assembly38.fasta"
WORKDIR = "/home/ozan/grad_project/cosap-nx/workdir_homework"


def main():
    bam = BamReader(filename=BAM_PATH, name="NA12878_exome")

    caller = VariantCaller(
        library="deepvariant",
        normal_sample=bam,
        params={"germline_sample_name": "NA12878_exome", "model_type": "WES"},
        gvcf=True,
    )

    pipeline = Pipeline(ref_fasta=REF_PATH)
    pipeline.add(caller)

    config_path = pipeline.build(workdir=WORKDIR, cpus=8, memory="14 GB")

    runner = PipelineRunner()
    output_dir = runner.run_pipeline(config_path)

    print(f"\nDone. Outputs:")
    print(f"  VCF:  {WORKDIR}/results/vcf/deepvariant/NA12878_exome.vcf.gz")
    print(f"  gVCF: {WORKDIR}/results/vcf/deepvariant/NA12878_exome.g.vcf.gz")


if __name__ == "__main__":
    main()
