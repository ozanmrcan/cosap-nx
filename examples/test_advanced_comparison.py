#!/usr/bin/env python3
"""
Test advanced comparison features: get_variants() and create_intersection_bed().

This demonstrates COSAP's additional utility methods for variant analysis.
"""

from cosap_nx import VariantComparator

# VCF outputs from previous pipeline runs
DEEPVARIANT_VCF = "./workdir/results/vcf/deepvariant/NA12878_chr20.vcf.gz"
HAPLOTYPECALLER_VCF = "./workdir_hc/results/vcf/haplotypecaller/NA12878_chr20.vcf.gz"

COMPARISON_DIR = "./comparison_results"


def main():
    print("=" * 70)
    print("Testing Advanced Comparison Features")
    print("=" * 70)

    # Create comparator
    comparator = VariantComparator(
        variant_files=[
            {"path": DEEPVARIANT_VCF, "caller": "deepvariant", "mapper": "bwa"},
            {"path": HAPLOTYPECALLER_VCF, "caller": "haplotypecaller", "mapper": "bwa"},
        ]
    )

    # Test 1: get_variants()
    print("\n[1/4] Testing get_variants()...")
    variants_df = comparator.get_variants()
    print(f"      Loaded {len(variants_df)} variant records")
    print(f"      Columns: {list(variants_df.columns)}")
    print(f"\n      First 5 variants:")
    print(variants_df.head())
    
    # Show counts per pipeline
    print(f"\n      Variants per pipeline:")
    print(variants_df['pipeline'].value_counts())

    # Test 2: create_intersection_bed() - Shared variants
    print("\n[2/4] Creating BED file for shared variants...")
    comparator.create_intersection_bed(
        "bwa_deepvariant & bwa_haplotypecaller",
        f"{COMPARISON_DIR}/shared_variants.bed"
    )
    print(f"      Saved: {COMPARISON_DIR}/shared_variants.bed")

    # Test 3: create_intersection_bed() - DeepVariant unique
    print("\n[3/4] Creating BED file for DeepVariant-only variants...")
    comparator.create_intersection_bed(
        "bwa_deepvariant ~ bwa_haplotypecaller",
        f"{COMPARISON_DIR}/deepvariant_unique.bed"
    )
    print(f"      Saved: {COMPARISON_DIR}/deepvariant_unique.bed")

    # Test 4: create_intersection_bed() - HaplotypeCaller unique
    print("\n[4/4] Creating BED file for HaplotypeCaller-only variants...")
    comparator.create_intersection_bed(
        "bwa_haplotypecaller ~ bwa_deepvariant",
        f"{COMPARISON_DIR}/haplotypecaller_unique.bed"
    )
    print(f"      Saved: {COMPARISON_DIR}/haplotypecaller_unique.bed")

    print("\n" + "=" * 70)
    print("Testing Complete!")
    print("=" * 70)

    # Show BED file examples
    print("\nBED File Preview (shared variants, first 5 lines):")
    import subprocess
    result = subprocess.run(
        ["head", "-5", f"{COMPARISON_DIR}/shared_variants.bed"],
        capture_output=True,
        text=True
    )
    print(result.stdout)

    print("\nUsage Examples:")
    print("  - Load in IGV for visual inspection")
    print("  - Use with bedtools for downstream analysis")
    print("  - Filter VCF files to specific variant sets")


if __name__ == "__main__":
    main()

