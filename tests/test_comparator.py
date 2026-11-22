import json

import pandas as pd
import pytest

from conftest import DEEPVARIANT_VCF, HAPLOTYPECALLER_VCF, TRUTH_VCF
from cosap_nx import PipelineRunner, VariantComparator

# Expected values from the 2026-11-21 baseline run.


def test_metrics_vs_truth(comparator, tmp_path):
    out = comparator.compute_metrics_vs_truth(truth_vcf=str(TRUTH_VCF), output_dir=str(tmp_path))
    metrics = json.loads(out.read_text())

    assert metrics == {
        "bwa_deepvariant": {
            "true_positives": 219,
            "false_positives": 69,
            "false_negatives": 2,
            "precision": 0.7604,
            "recall": 0.991,
            "f1": 0.8605,
        },
        "bwa_haplotypecaller": {
            "true_positives": 219,
            "false_positives": 28,
            "false_negatives": 2,
            "precision": 0.8866,
            "recall": 0.991,
            "f1": 0.9359,
        },
    }


def test_overlap(comparator, tmp_path):
    out = comparator.compute_overlap(output_dir=str(tmp_path))
    overlap = json.loads(out.read_text())

    assert overlap["bwa_deepvariant_vs_bwa_haplotypecaller"] == {
        "pipeline1": "bwa_deepvariant",
        "pipeline2": "bwa_haplotypecaller",
        "only_in_pipeline1": 48,
        "only_in_pipeline2": 7,
        "shared": 240,
        "jaccard_similarity": 0.8136,
    }


def test_statistics(comparator, tmp_path):
    out = comparator.compute_statistics(output_dir=str(tmp_path))
    counts = pd.read_csv(out, sep="\t").set_index("pipeline")

    # total = snps + indels, so HC's multiallelic SNP+indel site counts twice (247 records -> 248)
    assert counts.loc["bwa_deepvariant", ["total_variants", "snps", "indels"]].tolist() == [288, 216, 72]
    assert counts.loc["bwa_haplotypecaller", ["total_variants", "snps", "indels"]].tolist() == [248, 196, 52]


@pytest.mark.parametrize(
    "query, expected_lines",
    [
        ("bwa_deepvariant & bwa_haplotypecaller", 240),
        ("bwa_deepvariant ~ bwa_haplotypecaller", 48),
        ("bwa_haplotypecaller ~ bwa_deepvariant", 7),
    ],
)
def test_intersection_bed(comparator, tmp_path, query, expected_lines):
    bed = tmp_path / "out.bed"
    comparator.create_intersection_bed(query, str(bed))
    rows = [line.split("\t") for line in bed.read_text().splitlines()]

    assert len(rows) == expected_lines
    assert all(int(end) - int(start) == 80 for _, start, end in rows)


# variants.json is what the workflow writes into results/ after variant calling


@pytest.fixture
def results(tmp_path):
    entries = [
        {"name": "s_bwa_deepvariant", "path": str(DEEPVARIANT_VCF), "caller": "deepvariant", "mapper": "bwa", "sample": "s"},
        {"name": "s_bwa_haplotypecaller", "path": str(HAPLOTYPECALLER_VCF), "caller": "haplotypecaller", "mapper": "bwa", "sample": "s"},
        {"name": "t_bwa_deepvariant", "path": str(DEEPVARIANT_VCF), "caller": "deepvariant", "mapper": "bwa", "sample": "t"},
    ]
    (tmp_path / "variants.json").write_text(json.dumps(entries))
    return tmp_path


def test_from_results_one_sample(results):
    comparator = VariantComparator.from_results(str(results), sample="s")

    assert comparator.pipeline_names == ["bwa_deepvariant", "bwa_haplotypecaller"]


def test_from_results_same_label_twice(results):
    # s and t both have bwa_deepvariant
    with pytest.raises(ValueError):
        VariantComparator.from_results(str(results))


def test_from_results_without_manifest(tmp_path):
    with pytest.raises(FileNotFoundError):
        VariantComparator.from_results(str(tmp_path))


def test_compare_results_per_sample(results):
    out = PipelineRunner.compare_results(results, truth_vcf=str(TRUTH_VCF))

    assert sorted(p.name for p in out.iterdir()) == ["s", "t"]
    metrics = json.loads((out / "s" / "metrics_vs_truth.json").read_text())
    assert metrics["bwa_deepvariant"]["f1"] == 0.8605
    assert (out / "s" / "overlap_summary.json").exists()
    assert (out / "t" / "variant_counts.tsv").exists()
