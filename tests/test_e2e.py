import json
import shutil
import subprocess
from pathlib import Path

import pytest

from conftest import TRUTH_VCF
from cosap_nx import BamReader, FastqReader, Pipeline, PipelineRunner, VariantCaller, VariantComparator

# Full Nextflow runs on the chr20 inputs in ../data. Run with: pytest -m slow

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


def run(caller, workdir, cpus, memory):
    params = Pipeline(ref_fasta=str(REF)).add(caller).build(workdir=str(workdir), cpus=cpus, memory=memory)
    PipelineRunner().run_pipeline(params, resume=False)
    return workdir / "results" / "vcf" / caller.library / f"{SAMPLE}.vcf.gz"


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


def test_bam_deepvariant(tmp_path):
    caller = VariantCaller(
        library="deepvariant",
        normal_sample=BamReader(str(BAM), name=SAMPLE),
        params={"germline_sample_name": SAMPLE, "model_type": "WGS"},
    )
    vcf = run(caller, tmp_path / "work", cpus=4, memory="12 GB")

    assert count_records(vcf) == 288
    assert count_records(vcf, pass_only=True) == 238
    assert f1_vs_truth(vcf, "deepvariant", tmp_path) == 0.8605


def test_bam_haplotypecaller(tmp_path):
    caller = VariantCaller(
        library="haplotypecaller",
        normal_sample=BamReader(str(BAM), name=SAMPLE),
        params={"germline_sample_name": SAMPLE},
    )
    vcf = run(caller, tmp_path / "work", cpus=4, memory="8 GB")

    assert count_records(vcf) == 247
    assert f1_vs_truth(vcf, "haplotypecaller", tmp_path) == 0.9359


def test_fastq_deepvariant(tmp_path):
    reads = [
        FastqReader(str(FASTQ_R1), read=1, name=SAMPLE),
        FastqReader(str(FASTQ_R2), read=2, name=SAMPLE),
    ]
    caller = VariantCaller(library="deepvariant", normal_sample=reads, params={"model_type": "WGS"})
    vcf = run(caller, tmp_path / "work", cpus=8, memory="12 GB")

    assert (tmp_path / "work" / "results" / "bam" / f"{SAMPLE}.bam").exists()
    assert count_records(vcf) == 24
    assert count_records(vcf, pass_only=True) == 23
