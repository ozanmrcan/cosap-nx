import json
import shutil
import subprocess
from pathlib import Path

import pytest

from conftest import TRUTH_VCF
from cosap_nx import BamReader, FastqReader, Mapper, Pipeline, PipelineRunner, VariantCaller, VariantComparator

# One full Nextflow run on the chr20 inputs in ../data. Run with: pytest -m slow

INPUTS = Path(__file__).parent.parent.parent / "data"
REF = INPUTS / "refs" / "ucsc.hg19.chr20.unittest.fasta"
BAM = INPUTS / "bam" / "NA12878_S1.chr20.10_10p1mb.bam"
FASTQ_R1 = INPUTS / "fastq" / "NA12878_chr20_R1.fastq.gz"
FASTQ_R2 = INPUTS / "fastq" / "NA12878_chr20_R2.fastq.gz"
SAMPLE = "NA12878_chr20"

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not REF.exists(), reason=f"chr20 inputs not found in {INPUTS}"),
    pytest.mark.skipif(not shutil.which("docker"), reason="docker not found"),
    pytest.mark.skipif(not shutil.which("nextflow"), reason="nextflow not found"),
]


@pytest.fixture(scope="module")
def results(tmp_path_factory):
    """BAM -> DeepVariant/HaplotypeCaller and FASTQ -> BWA/BWA2 -> DeepVariant/HaplotypeCaller in one run."""
    bam = BamReader(str(BAM), name=f"{SAMPLE}_bam")
    reads = [
        FastqReader(str(FASTQ_R1), read=1, name=SAMPLE),
        FastqReader(str(FASTQ_R2), read=2, name=SAMPLE),
    ]
    bwa = Mapper(library="bwa", input_step=reads, params={"read_groups": {"SM": SAMPLE}})
    bwa2 = Mapper(library="bwa2", input_step=reads, params={"read_groups": {"SM": SAMPLE}})

    pipeline = (
        Pipeline(ref_fasta=str(REF))
        .add(VariantCaller(library="deepvariant", germline=bam, params={"model_type": "WGS"}))
        .add(VariantCaller(library="haplotypecaller", germline=bam))
        .add(VariantCaller(library="deepvariant", germline=bwa, params={"model_type": "WGS"}))
        .add(VariantCaller(library="haplotypecaller", germline=bwa))
        .add(VariantCaller(library="deepvariant", germline=bwa2, params={"model_type": "WGS"}))
        .add(VariantCaller(library="haplotypecaller", germline=bwa2))
    )
    workdir = tmp_path_factory.mktemp("work")
    params = pipeline.build(workdir=str(workdir), cpus=8, memory="12 GB")
    return PipelineRunner().run_pipeline(params, resume=False)


def count_records(vcf, pass_only=False):
    cmd = ["bcftools", "view", "-H", str(vcf)]
    if pass_only:
        cmd[2:2] = ["-f", "PASS"]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    return len(out.splitlines())


def f1_vs_truth(vcf, caller, tmp_path):
    comparator = VariantComparator(variant_files=[{"path": str(vcf), "caller": caller, "mapper": "bwa"}])
    out = comparator.compute_metrics_vs_truth(truth_vcf=str(TRUTH_VCF), output_dir=str(tmp_path / "cmp"))
    return json.loads(out.read_text())[f"bwa_{caller}"]["f1"]


def test_bam_deepvariant(results, tmp_path):
    vcf = results / "vcf" / "deepvariant" / f"{SAMPLE}_bam_deepvariant.vcf.gz"

    assert count_records(vcf) == 288
    assert count_records(vcf, pass_only=True) == 238
    assert f1_vs_truth(vcf, "deepvariant", tmp_path) == 0.8605


def test_bam_haplotypecaller(results, tmp_path):
    vcf = results / "vcf" / "haplotypecaller" / f"{SAMPLE}_bam_haplotypecaller.vcf.gz"

    assert count_records(vcf) == 247
    assert f1_vs_truth(vcf, "haplotypecaller", tmp_path) == 0.9359


def test_fastq_deepvariant(results):
    vcf = results / "vcf" / "deepvariant" / f"{SAMPLE}_bwa_deepvariant.vcf.gz"

    assert (results / "bam" / "bwa" / f"{SAMPLE}_bwa.bam").exists()
    assert count_records(vcf) == 24
    assert count_records(vcf, pass_only=True) == 23


def test_fastq_haplotypecaller(results):
    vcf = results / "vcf" / "haplotypecaller" / f"{SAMPLE}_bwa_haplotypecaller.vcf.gz"

    assert count_records(vcf) == 16


def test_fastq_bwa2_deepvariant(results):
    vcf = results / "vcf" / "deepvariant" / f"{SAMPLE}_bwa2_deepvariant.vcf.gz"

    assert (results / "bam" / "bwa2" / f"{SAMPLE}_bwa2.bam").exists()
    assert count_records(vcf) == 24


def test_fastq_bwa2_haplotypecaller(results):
    vcf = results / "vcf" / "haplotypecaller" / f"{SAMPLE}_bwa2_haplotypecaller.vcf.gz"

    assert count_records(vcf) == 16
