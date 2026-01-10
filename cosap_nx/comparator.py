"""Variant comparison functionality matching COSAP API."""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Optional
import subprocess
import json
import pandas as pd


@dataclass
class VariantFile:
    """Represents a variant file with metadata."""
    path: str
    caller: str
    mapper: str


class VariantComparator:
    """
    Compare multiple VCF files from different pipelines.

    API matches COSAP's VariantComparator exactly.

    Example:
        comparator = VariantComparator(
            variant_files=[
                {"path": "vcf1.vcf.gz", "caller": "deepvariant", "mapper": "bwa"},
                {"path": "vcf2.vcf.gz", "caller": "haplotypecaller", "mapper": "bwa"},
            ]
        )
        comparator.compute_statistics(output_dir="./comparison_results")
        comparator.compute_overlap(output_dir="./comparison_results")
        comparator.compute_metrics_vs_truth(
            truth_vcf="truth.vcf.gz",
            output_dir="./comparison_results"
        )
    """

    def __init__(self, variant_files: List[Dict[str, str]], bed_file: Optional[str] = None):
        """
        Args:
            variant_files: List of dicts with keys:
                - path: Path to VCF file
                - caller: Caller name (e.g., "deepvariant", "haplotypecaller")
                - mapper: Mapper name (e.g., "bwa")
            bed_file: Optional BED file for region filtering (not implemented in v0.2.0)
        """
        self.variant_files = [VariantFile(**vf) for vf in variant_files]
        self.bed_file = bed_file
        self.pipeline_names = [f"{vf.mapper}_{vf.caller}" for vf in self.variant_files]

        # Validate VCF files exist
        for vf in self.variant_files:
            if not Path(vf.path).exists():
                raise FileNotFoundError(f"VCF file not found: {vf.path}")

    def compute_statistics(self, output_dir: str) -> Path:
        """
        Compute variant statistics for each VCF.

        Uses bcftools stats to count SNPs, indels, and total variants.

        Args:
            output_dir: Directory to save output files

        Returns:
            Path to variant_counts.tsv
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        stats = []
        for vf, pipeline_name in zip(self.variant_files, self.pipeline_names):
            # Use bcftools stats to count variants
            cmd = ["bcftools", "stats", vf.path]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            # Parse bcftools stats output
            snp_count = 0
            indel_count = 0
            for line in result.stdout.split('\n'):
                if line.startswith('SN') and 'number of SNPs:' in line:
                    snp_count = int(line.split('\t')[3])
                elif line.startswith('SN') and 'number of indels:' in line:
                    indel_count = int(line.split('\t')[3])

            stats.append({
                'pipeline': pipeline_name,
                'caller': vf.caller,
                'mapper': vf.mapper,
                'total_variants': snp_count + indel_count,
                'snps': snp_count,
                'indels': indel_count
            })

        # Write to TSV
        df = pd.DataFrame(stats)
        output_file = output_path / "variant_counts.tsv"
        df.to_csv(output_file, sep='\t', index=False)

        return output_file

    def compute_overlap(self, output_dir: str) -> Path:
        """
        Compute pairwise overlap between VCFs using bcftools isec.

        Calculates:
        - Variants unique to each pipeline
        - Shared variants
        - Jaccard similarity coefficient

        Args:
            output_dir: Directory to save output files

        Returns:
            Path to overlap_summary.json
        """
        output_path = Path(output_dir)
        isec_dir = output_path / "isec_outputs"
        isec_dir.mkdir(parents=True, exist_ok=True)

        overlap_results = {}

        # Pairwise comparison
        for i in range(len(self.variant_files)):
            for j in range(i + 1, len(self.variant_files)):
                vf1 = self.variant_files[i]
                vf2 = self.variant_files[j]
                pipeline1 = self.pipeline_names[i]
                pipeline2 = self.pipeline_names[j]

                pair_dir = isec_dir / f"{pipeline1}_vs_{pipeline2}"
                pair_dir.mkdir(exist_ok=True)

                # Run bcftools isec
                # Output files: 0000.vcf (only in vcf1), 0001.vcf (only in vcf2),
                #               0002.vcf (shared, from vcf1), 0003.vcf (shared, from vcf2)
                cmd = ["bcftools", "isec", vf1.path, vf2.path, "-p", str(pair_dir)]
                subprocess.run(cmd, check=True, capture_output=True)

                # Count variants in each output file
                def count_vcf_variants(vcf_path):
                    if not vcf_path.exists():
                        return 0
                    cmd = ["bcftools", "view", "-H", str(vcf_path)]
                    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                    return len([line for line in result.stdout.split('\n') if line.strip()])

                only_in_1 = count_vcf_variants(pair_dir / "0000.vcf")
                only_in_2 = count_vcf_variants(pair_dir / "0001.vcf")
                shared = count_vcf_variants(pair_dir / "0002.vcf")

                # Calculate Jaccard similarity: |A ∩ B| / |A ∪ B|
                total_union = only_in_1 + only_in_2 + shared
                jaccard = shared / total_union if total_union > 0 else 0

                overlap_results[f"{pipeline1}_vs_{pipeline2}"] = {
                    "pipeline1": pipeline1,
                    "pipeline2": pipeline2,
                    "only_in_pipeline1": only_in_1,
                    "only_in_pipeline2": only_in_2,
                    "shared": shared,
                    "jaccard_similarity": round(jaccard, 4)
                }

        # Write summary
        output_file = output_path / "overlap_summary.json"
        with open(output_file, 'w') as f:
            json.dump(overlap_results, f, indent=2)

        return output_file

    def compute_metrics_vs_truth(self, truth_vcf: str, output_dir: str) -> Path:
        """
        Compute precision/recall/F1 vs ground truth VCF.

        Uses bcftools isec to identify:
        - True Positives (TP): Variants in both query and truth
        - False Positives (FP): Variants in query but not in truth
        - False Negatives (FN): Variants in truth but not in query

        Metrics:
        - Precision = TP / (TP + FP)
        - Recall = TP / (TP + FN)
        - F1 = 2 * (Precision * Recall) / (Precision + Recall)

        Note: Automatically handles chromosome naming differences (chr20 vs 20)
        by attempting to normalize the truth VCF if needed.

        Args:
            truth_vcf: Path to ground truth VCF file
            output_dir: Directory to save output files

        Returns:
            Path to metrics_vs_truth.json
        """
        if not Path(truth_vcf).exists():
            raise FileNotFoundError(f"Truth VCF not found: {truth_vcf}")

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Check if chromosome naming matches between query and truth
        normalized_truth = self._normalize_chromosome_names(truth_vcf, output_path)

        metrics = {}

        for vf, pipeline_name in zip(self.variant_files, self.pipeline_names):
            # Use bcftools isec to find TP, FP, FN
            temp_dir = output_path / "temp_isec" / pipeline_name
            temp_dir.mkdir(parents=True, exist_ok=True)

            # Compare query VCF vs truth VCF (using normalized truth if needed)
            cmd = ["bcftools", "isec", vf.path, normalized_truth, "-p", str(temp_dir)]
            subprocess.run(cmd, check=True, capture_output=True)

            def count_vcf_variants(vcf_path):
                if not vcf_path.exists():
                    return 0
                cmd = ["bcftools", "view", "-H", str(vcf_path)]
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                return len([line for line in result.stdout.split('\n') if line.strip()])

            fp = count_vcf_variants(temp_dir / "0000.vcf")  # Only in query (FP)
            fn = count_vcf_variants(temp_dir / "0001.vcf")  # Only in truth (FN)
            tp = count_vcf_variants(temp_dir / "0002.vcf")  # Shared (TP)

            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

            metrics[pipeline_name] = {
                "true_positives": tp,
                "false_positives": fp,
                "false_negatives": fn,
                "precision": round(precision, 4),
                "recall": round(recall, 4),
                "f1": round(f1, 4)
            }

        # Write metrics
        output_file = output_path / "metrics_vs_truth.json"
        with open(output_file, 'w') as f:
            json.dump(metrics, f, indent=2)

        return output_file

    def _normalize_chromosome_names(self, truth_vcf: str, output_dir: Path) -> str:
        """
        Normalize chromosome names in truth VCF to match query VCFs.

        Checks if the first query VCF uses 'chr' prefix. If yes and truth doesn't,
        creates a normalized version of the truth VCF with 'chr' prefix added.

        Args:
            truth_vcf: Path to original truth VCF
            output_dir: Directory to save normalized VCF

        Returns:
            Path to normalized truth VCF (or original if no normalization needed)
        """
        # Check chromosome naming in first query VCF
        first_vcf = self.variant_files[0].path
        cmd = ["bcftools", "view", "-H", first_vcf]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        first_line = result.stdout.split('\n')[0] if result.stdout else ""
        query_has_chr = first_line.startswith('chr') if first_line else False

        # Check chromosome naming in truth VCF
        cmd = ["bcftools", "view", "-H", truth_vcf]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        first_line = result.stdout.split('\n')[0] if result.stdout else ""
        truth_has_chr = first_line.startswith('chr') if first_line else False

        # If naming matches, return original
        if query_has_chr == truth_has_chr:
            return truth_vcf

        # Need to normalize - add 'chr' prefix to truth VCF
        normalized_dir = output_dir / "normalized_truth"
        normalized_dir.mkdir(exist_ok=True)
        normalized_vcf = normalized_dir / "truth_normalized.vcf.gz"

        if query_has_chr and not truth_has_chr:
            # Add 'chr' prefix to truth VCF chromosomes
            cmd = f"bcftools annotate --rename-chrs <(echo '20 chr20') {truth_vcf} -Oz -o {normalized_vcf}"
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
            # Index the normalized VCF
            subprocess.run(["bcftools", "index", "-t", str(normalized_vcf)], check=True)
            return str(normalized_vcf)
        elif not query_has_chr and truth_has_chr:
            # Remove 'chr' prefix from truth VCF chromosomes
            cmd = f"bcftools annotate --rename-chrs <(echo 'chr20 20') {truth_vcf} -Oz -o {normalized_vcf}"
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
            # Index the normalized VCF
            subprocess.run(["bcftools", "index", "-t", str(normalized_vcf)], check=True)
            return str(normalized_vcf)

        return truth_vcf

