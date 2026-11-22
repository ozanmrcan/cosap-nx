import json
import shutil
import subprocess
from pathlib import Path

import pytest

from cosap_nx import MDUP, FastqReader, Mapper, Pipeline, Recalibrator, Trimmer, VariantCaller

# Runs nf/main.nf with -stub-run: no tools or containers, only checks that
# every step runs and gets the right inputs. Run with: pytest -m stub

MAIN_NF = Path(__file__).parent.parent / "nf" / "main.nf"

pytestmark = [
    pytest.mark.stub,
    pytest.mark.skipif(not shutil.which("nextflow"), reason="nextflow not found"),
]


def stub_run(pipeline, tmp_path):
    (tmp_path / "ref.fasta").touch()
    (tmp_path / "ref.fasta.fai").touch()
    params = pipeline.set_reference(str(tmp_path / "ref.fasta")).build(workdir=str(tmp_path / "work"), cpus=1, memory="1 GB")
    cmd = ["nextflow", "run", str(MAIN_NF), "-params-file", params, "-stub-run"]
    subprocess.run(cmd, cwd=tmp_path / "work", check=True, capture_output=True)
    return tmp_path / "work" / "results"


def reads(tmp_path, sample):
    for i in [1, 2]:
        (tmp_path / f"{sample}_R{i}.fastq.gz").touch()
    return [FastqReader(str(tmp_path / f"{sample}_R{i}.fastq.gz"), read=i) for i in [1, 2]]


def test_matrix(tmp_path):
    pipeline = Pipeline()
    for m in ["bwa", "bwa2", "bowtie2"]:
        mapped = Mapper(library=m, input_step=reads(tmp_path, "normal"), params={"read_groups": {"SM": "normal"}})
        for c in ["deepvariant", "haplotypecaller", "strelka"]:
            pipeline.add(VariantCaller(library=c, germline=mapped))
    results = stub_run(pipeline, tmp_path)

    vcfs = sorted(p.name for p in results.rglob("*.vcf*") if not p.name.endswith(".tbi"))
    assert vcfs == [
        "normal_bowtie2_deepvariant.vcf.gz",
        "normal_bowtie2_haplotypecaller.vcf.gz",
        "normal_bowtie2_strelka.vcf",
        "normal_bwa2_deepvariant.vcf.gz",
        "normal_bwa2_haplotypecaller.vcf.gz",
        "normal_bwa2_strelka.vcf",
        "normal_bwa_deepvariant.vcf.gz",
        "normal_bwa_haplotypecaller.vcf.gz",
        "normal_bwa_strelka.vcf",
    ]

    # variants.json lists every call with its labels, at its published path
    entries = json.loads((results / "variants.json").read_text())
    assert sorted((e["mapper"], e["caller"]) for e in entries) == sorted(
        (m, c) for m in ["bwa", "bwa2", "bowtie2"] for c in ["deepvariant", "haplotypecaller", "strelka"]
    )
    assert all(Path(e["path"]).exists() and e["sample"] == "normal" for e in entries)


def test_somatic_preprocessing_chain(tmp_path):
    def chain(sample):
        mapped = Mapper(library="bwa", input_step=Trimmer(input_step=reads(tmp_path, sample)),
                        params={"read_groups": {"SM": sample}})
        return Recalibrator(input_step=MDUP(input_step=mapped))

    pipeline = Pipeline().add(VariantCaller(library="mutect", germline=chain("normal"), tumor=chain("tumor")))
    results = stub_run(pipeline, tmp_path)

    # The placeholder writes "<tool> on: <input files>", normal first
    mutect = results / "variant_calling" / "normal_bwa_tumor_bwa_mutect.vcf"
    assert mutect.read_text().strip() == "mutect on: normal_bwa.bam tumor_bwa.bam"
