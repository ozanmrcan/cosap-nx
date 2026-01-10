#!/usr/bin/env python3
"""
Compare DeepVariant vs HaplotypeCaller on chr20 test data.

This example demonstrates COSAP-NX's comparison workflow:
1. Run DeepVariant pipeline (generates VCF)
2. Run HaplotypeCaller pipeline (generates VCF)
3. Compare the two VCFs using VariantComparator

Note: This script assumes you've already run both pipelines.
To run the pipelines first:
    python examples/germline_deepvariant.py
    python examples/germline_haplotypecaller.py
"""

from cosap_nx import VariantComparator

# VCF outputs from previous pipeline runs
DEEPVARIANT_VCF = "./workdir/results/vcf/deepvariant/NA12878_chr20.vcf.gz"
HAPLOTYPECALLER_VCF = "./workdir_hc/results/vcf/haplotypecaller/NA12878_chr20.vcf.gz"

# Optional: truth VCF for precision/recall
TRUTH_VCF = "/home/ozan/grad_project/data/benchmark/truth_chr20_10-10.1mb.vcf.gz"

COMPARISON_DIR = "./comparison_results"


def main():
    print("=" * 70)
    print("COSAP-NX Variant Comparison: DeepVariant vs HaplotypeCaller")
    print("=" * 70)

    # Create comparator with COSAP-compatible API
    comparator = VariantComparator(
        variant_files=[
            {"path": DEEPVARIANT_VCF, "caller": "deepvariant", "mapper": "bwa"},
            {"path": HAPLOTYPECALLER_VCF, "caller": "haplotypecaller", "mapper": "bwa"},
        ]
    )

    print("\n[1/3] Computing variant statistics...")
    stats_file = comparator.compute_statistics(output_dir=COMPARISON_DIR)
    print(f"      Saved: {stats_file}")

    print("\n[2/3] Computing pairwise overlap...")
    overlap_file = comparator.compute_overlap(output_dir=COMPARISON_DIR)
    print(f"      Saved: {overlap_file}")

    print("\n[3/3] Computing metrics vs truth VCF...")
    metrics_file = comparator.compute_metrics_vs_truth(
        truth_vcf=TRUTH_VCF,
        output_dir=COMPARISON_DIR
    )
    print(f"      Saved: {metrics_file}")

    print("\n" + "=" * 70)
    print(f"Done! All results saved to: {COMPARISON_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()

