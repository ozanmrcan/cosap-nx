import json

import pytest

from cosap_nx import MDUP, BamReader, FastqReader, Mapper, Pipeline, Recalibrator, Trimmer, VariantCaller

# build() only checks that input files exist, so empty placeholder files are enough.


@pytest.fixture
def files(tmp_path):
    names = ["ref.fasta", "ref.fasta.fai", "sample.bam", "sample.bam.bai", "sample_R1.fastq.gz", "sample_R2.fastq.gz"]
    for name in names:
        (tmp_path / name).touch()
    return {name: str(tmp_path / name) for name in names}


def fastqs(files):
    return [
        FastqReader(files["sample_R1.fastq.gz"], read=1),
        FastqReader(files["sample_R2.fastq.gz"], read=2),
    ]


def mapper(files, library="bwa", read_groups={"SM": "sample"}):
    return Mapper(library=library, input_step=fastqs(files), params={"read_groups": read_groups})


def build(files, tmp_path, *steps):
    pipeline = Pipeline(ref_fasta=files["ref.fasta"])
    for step in steps:
        pipeline.add(step)
    path = pipeline.build(workdir=str(tmp_path / "work"))
    assert (tmp_path / "work" / "params.json").samefile(path)
    return json.load(open(path))


# The JSON format is the contract with nf/main.nf, so its shape is checked here.


def test_bam_caller(files, tmp_path):
    caller = VariantCaller(library="deepvariant", germline=BamReader(files["sample.bam"]))
    config = build(files, tmp_path, caller)

    assert config["inputs"] == {"sample": {"type": "bam", "bam": files["sample.bam"]}}
    assert config["variant_calling"]["sample_deepvariant"]["normal"] == "inputs/sample"


def test_mapper_feeds_callers(files, tmp_path):
    bwa = mapper(files)
    dv = VariantCaller(library="deepvariant", germline=bwa)
    hc = VariantCaller(library="haplotypecaller", germline=bwa)
    config = build(files, tmp_path, bwa, dv, hc)

    assert config["mapping"]["sample_bwa"]["input"] == "inputs/sample"
    assert config["variant_calling"]["sample_bwa_deepvariant"]["normal"] == "mapping/sample_bwa"
    assert config["variant_calling"]["sample_bwa_haplotypecaller"]["normal"] == "mapping/sample_bwa"


def test_upstream_steps_are_added(files, tmp_path):
    # Only the caller is added, the mapper comes with it
    caller = VariantCaller(library="deepvariant", germline=mapper(files))
    config = build(files, tmp_path, caller)

    assert "sample_bwa" in config["mapping"]


def test_matrix(files, tmp_path):
    steps = []
    for m in ["bwa", "bowtie2"]:
        mapped = mapper(files, library=m)
        steps += [VariantCaller(library=c, germline=mapped) for c in ["deepvariant", "haplotypecaller", "strelka"]]
    config = build(files, tmp_path, *steps)

    assert len(config["inputs"]) == 1
    assert len(config["mapping"]) == 2
    assert len(config["variant_calling"]) == 6


def test_somatic_preprocessing_chain(files, tmp_path):
    # Same shape as COSAP's examples/example_run.py
    def chain(sample):
        reads = [FastqReader(files["sample_R1.fastq.gz"], read=1, name=sample),
                 FastqReader(files["sample_R2.fastq.gz"], read=2, name=sample)]
        mapped = Mapper(library="bwa", input_step=Trimmer(input_step=reads), params={"read_groups": {"SM": sample}})
        return Recalibrator(input_step=MDUP(input_step=mapped))

    mutect = VariantCaller(library="mutect", germline=chain("normal"), tumor=chain("tumor"))
    config = build(files, tmp_path, mutect)

    assert config["mapping"]["normal_bwa"]["input"] == "trim/normal"
    assert config["mdup"]["normal_bwa"] == {"input": "mapping/normal_bwa", "spark": True}
    assert config["calibrate"]["tumor_bwa"]["input"] == "mdup/tumor_bwa"
    assert config["variant_calling"]["normal_bwa_tumor_bwa_mutect"]["normal"] == "calibrate/normal_bwa"
    assert config["variant_calling"]["normal_bwa_tumor_bwa_mutect"]["tumor"] == "calibrate/tumor_bwa"


def test_resources_are_defaults(files, tmp_path):
    caller = VariantCaller(library="deepvariant", germline=BamReader(files["sample.bam"]))
    path = Pipeline(ref_fasta=files["ref.fasta"]).add(caller).build(workdir=str(tmp_path / "work"), cpus=2, memory="4 GB")

    assert json.load(open(path))["defaults"] == {"cpus": 2, "memory": "4 GB"}


def test_same_name_twice(files, tmp_path):
    a = VariantCaller(library="deepvariant", germline=BamReader(files["sample.bam"]), name="calls")
    b = VariantCaller(library="haplotypecaller", germline=BamReader(files["sample.bam"]), name="calls")
    with pytest.raises(ValueError):
        build(files, tmp_path, a, b)


def test_caller_needs_input():
    with pytest.raises(ValueError):
        VariantCaller(library="deepvariant")


def test_mapper_needs_sample_name(files):
    with pytest.raises(ValueError):
        mapper(files, read_groups={"ID": "run1"})


def test_fastq_reads_must_be_1_and_2(files):
    reads = [
        FastqReader(files["sample_R1.fastq.gz"], read=1),
        FastqReader(files["sample_R2.fastq.gz"], read=1),
    ]
    with pytest.raises(ValueError):
        Mapper(library="bwa", input_step=reads, params={"read_groups": {"SM": "s"}})


def test_missing_fastq(tmp_path):
    with pytest.raises(FileNotFoundError):
        FastqReader(str(tmp_path / "nope_R1.fastq.gz"), read=1)


def test_missing_bam_index(files, tmp_path):
    (tmp_path / "sample.bam.bai").unlink()
    caller = VariantCaller(library="deepvariant", germline=BamReader(files["sample.bam"]))
    with pytest.raises(FileNotFoundError):
        build(files, tmp_path, caller)


def test_missing_reference_index(files, tmp_path):
    (tmp_path / "ref.fasta.fai").unlink()
    caller = VariantCaller(library="deepvariant", germline=BamReader(files["sample.bam"]))
    with pytest.raises(FileNotFoundError):
        build(files, tmp_path, caller)
