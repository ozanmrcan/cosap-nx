"""Variant comparison functionality matching COSAP API."""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Optional
import subprocess
import json
import pandas as pd

# Visualization imports
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn3
import upsetplot


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

    def _load_variants_dataframe(self) -> pd.DataFrame:
        """
        Load variants from all VCF files into a pandas DataFrame.

        Uses bcftools query to efficiently extract variant positions without
        parsing full VCF files.

        Returns:
            DataFrame with columns: variant_id, pipeline
            variant_id format: CHROM-POS-REF-ALT
        """
        variants_data = []

        for vf, pipeline_name in zip(self.variant_files, self.pipeline_names):
            # Use bcftools query to extract variant IDs
            cmd = ["bcftools", "query", "-f", "%CHROM-%POS-%REF-%ALT\\n", vf.path]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            # Parse output into list of variant IDs
            variant_ids = [line.strip() for line in result.stdout.split('\n') if line.strip()]

            # Add to data list
            for variant_id in variant_ids:
                variants_data.append({
                    'variant_id': variant_id,
                    'pipeline': pipeline_name
                })

        return pd.DataFrame(variants_data)

    def draw_venn_diagram(self, output_file: str = None):
        """
        Draw Venn diagram automatically detecting 2 or 3 pipelines.

        Args:
            output_file: Path to save the plot (PNG). If None, uses 'venn_diagram.png'
        """
        if len(self.pipeline_names) == 2:
            self.draw_venn2_plot(
                self.pipeline_names[0],
                self.pipeline_names[1],
                output_file=output_file
            )
        elif len(self.pipeline_names) == 3:
            self.draw_venn3_plot(
                self.pipeline_names[0],
                self.pipeline_names[1],
                self.pipeline_names[2],
                output_file=output_file
            )
        else:
            raise ValueError(f"Venn diagrams support 2-3 pipelines, got {len(self.pipeline_names)}")

    def draw_venn2_plot(self, pipeline1: str, pipeline2: str, output_file: str = None):
        """
        Draw 2-way Venn diagram comparing two pipelines.

        Args:
            pipeline1: Name of first pipeline
            pipeline2: Name of second pipeline
            output_file: Path to save the plot (PNG). If None, uses 'venn2_diagram.png'
        """
        if output_file is None:
            output_file = "venn2_diagram.png"

        # Load variant data
        variants_df = self._load_variants_dataframe()

        # Get variant sets for each pipeline - exactly like COSAP
        venn2(
            subsets=(
                set(variants_df[variants_df['pipeline'] == pipeline1]['variant_id'].tolist()),
                set(variants_df[variants_df['pipeline'] == pipeline2]['variant_id'].tolist()),
            ),
            set_labels=(pipeline1, pipeline2)
        )
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

    def draw_venn3_plot(self, pipeline1: str, pipeline2: str, pipeline3: str, output_file: str = None):
        """
        Draw 3-way Venn diagram comparing three pipelines.

        Args:
            pipeline1: Name of first pipeline
            pipeline2: Name of second pipeline
            pipeline3: Name of third pipeline
            output_file: Path to save the plot (PNG). If None, uses 'venn3_diagram.png'
        """
        if output_file is None:
            output_file = "venn3_diagram.png"

        # Load variant data
        variants_df = self._load_variants_dataframe()

        # Get variant sets for each pipeline - exactly like COSAP
        venn3(
            subsets=(
                set(variants_df[variants_df['pipeline'] == pipeline1]['variant_id'].tolist()),
                set(variants_df[variants_df['pipeline'] == pipeline2]['variant_id'].tolist()),
                set(variants_df[variants_df['pipeline'] == pipeline3]['variant_id'].tolist()),
            ),
            set_labels=(pipeline1, pipeline2, pipeline3)
        )
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

    def draw_upset_plot(self, output_file: str = None):
        """
        Draw UpSet plot showing all set intersections.

        UpSet plots are better than Venn diagrams for visualizing intersections
        of more than 3 sets, but work well for 2-3 sets too.

        Args:
            output_file: Path to save the plot (PNG). If None, uses 'upset_plot.png'
        """
        if output_file is None:
            output_file = "upset_plot.png"

        # Load variant data - COSAP approach
        variants_df = self._load_variants_dataframe()

        # Following COSAP's exact approach (lines 61-105 in comparator.py)
        upset_df = variants_df.copy()
        upset_df["value"] = True
        upset_df = upset_df.pivot(
            index="variant_id", columns="pipeline", values="value"
        ).replace(pd.NA, False)

        # Group by all pipeline columns to get intersection counts
        upset_df = (
            upset_df.groupby(
                list(upset_df.columns)
            )
            .size()
            .reset_index(name='count')
        )

        # Set pipelines as index
        upset_df = upset_df.set_index(list(upset_df.columns[:-1]))

        # Create UpSet plot - exactly like COSAP (no customization)
        upset = upsetplot.UpSet(
            upset_df['count'],
            sort_by="-degree",
        )

        upset.plot()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

    def draw_similarity_plot(self, output_file: str = None):
        """
        Draw heatmap showing Jaccard similarity between all pipelines.

        Jaccard similarity = |A ∩ B| / |A ∪ B|
        Perfect similarity = 1.0, no overlap = 0.0

        Args:
            output_file: Path to save the plot (PNG). If None, uses 'similarity_heatmap.png'
        """
        if output_file is None:
            output_file = "similarity_heatmap.png"

        # Load variant data
        variants_df = self._load_variants_dataframe()

        # Following COSAP's exact approach (lines 107-117)
        sim_df = variants_df.copy()

        sim_df = sim_df.pivot_table(
            index="variant_id", values="pipeline", aggfunc=lambda x: ",".join(x)
        )["pipeline"].str.get_dummies(sep=",")

        # Compute Jaccard similarity matrix
        from sklearn.metrics import pairwise_distances
        jac_sim = 1 - pairwise_distances(sim_df.T, metric="hamming")
        jac_sim = pd.DataFrame(jac_sim, columns=sim_df.columns, index=sim_df.columns)

        # Create heatmap - exactly like COSAP
        sns.heatmap(jac_sim, annot=True, cmap="YlGnBu")
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

    def draw_precision_recall_plot(self, truth_vcf: str, output_file: str = None):
        """
        Draw scatter plot comparing precision vs recall for each pipeline.

        Requires that compute_metrics_vs_truth() has been called first,
        or will compute metrics on the fly.

        Args:
            truth_vcf: Path to ground truth VCF file
            output_file: Path to save the plot (PNG). If None, uses 'precision_recall.png'
        """
        if output_file is None:
            output_file = "precision_recall.png"

        # Compute metrics if not already done (using temp directory)
        import tempfile
        with tempfile.TemporaryDirectory() as temp_dir:
            metrics_file = self.compute_metrics_vs_truth(truth_vcf, temp_dir)
            with open(metrics_file, 'r') as f:
                metrics_data = json.load(f)

        # Extract precision and recall for each pipeline
        precision_recall_values = {}
        for pipeline_name, metrics in metrics_data.items():
            precision_recall_values[pipeline_name] = (metrics['precision'], metrics['recall'])

        # Create DataFrame - following COSAP (lines 163-165)
        precision_recall_df = pd.DataFrame(precision_recall_values).T
        precision_recall_df.columns = ["precision", "recall"]

        # Create scatter plot - exactly like COSAP (lines 167-173)
        g = sns.scatterplot(
            data=precision_recall_df,
            x="precision",
            y="recall",
            hue=precision_recall_df.index,
            s=150,
        )

        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

    def get_variants(self) -> pd.DataFrame:
        """
        Get the variants DataFrame for custom analysis.

        Returns:
            DataFrame with columns: variant_id, pipeline
        """
        return self._load_variants_dataframe()

    def create_intersection_bed(self, query: str, output_file: str):
        """
        Create a BED file with variant positions from set operations.

        Supports set operations on pipelines:
        - & (intersection): variants in both
        - | (union): variants in either
        - ~ (difference): variants in first but not second

        The BED file includes ±40bp around each variant position.

        Args:
            query: Set operation query (e.g., "bwa_deepvariant & bwa_haplotypecaller")
            output_file: Output BED file path

        Example:
            # Get variants unique to DeepVariant
            comparator.create_intersection_bed(
                "bwa_deepvariant ~ bwa_haplotypecaller",
                "deepvariant_unique.bed"
            )
        """
        # Load variant data
        variants_df = self._load_variants_dataframe()

        # Set operation definitions (following COSAP's approach)
        operations = {
            "|": {"precedence": 1, "function": lambda x, y: x.union(y)},
            "&": {"precedence": 2, "function": lambda x, y: x.intersection(y)},
            "~": {"precedence": 3, "function": lambda x, y: x.difference(y)},
        }

        def shunting_yard(query):
            """Convert infix to postfix notation using Shunting Yard algorithm."""
            output = []
            stack = []
            for token in query.split():
                if token in operations:
                    while (
                        stack
                        and stack[-1] in operations
                        and operations[token]["precedence"]
                        <= operations[stack[-1]]["precedence"]
                    ):
                        output.append(stack.pop())
                    stack.append(token)
                else:
                    output.append(token)
            while stack:
                output.append(stack.pop())
            return output

        postfix = shunting_yard(query)

        # Evaluate postfix expression
        stack = []
        for token in postfix:
            if token in operations:
                y = stack.pop()
                x = stack.pop()

                # Convert pipeline names to variant sets
                x = (
                    set(variants_df[variants_df["pipeline"] == x]["variant_id"].tolist())
                    if isinstance(x, str)
                    else x
                )
                y = (
                    set(variants_df[variants_df["pipeline"] == y]["variant_id"].tolist())
                    if isinstance(y, str)
                    else y
                )

                results = operations[token]["function"](x, y)
                stack.append(results)
            else:
                stack.append(token)

        results = stack[0]

        # Create BED file from variant IDs
        bed_data = []
        for variant in results:
            chrom, pos, ref, alt = variant.split("-")
            start = int(pos) - 40  # ±40bp around variant
            end = int(pos) + 40
            bed_data.append({"CHROM": chrom, "START": start, "END": end})

        bed_df = pd.DataFrame(bed_data)
        bed_df.to_csv(output_file, sep="\t", index=False, header=False)





