import pytest

from cosap_nx import BamReader, FastqReader, Mapper, Pipeline, VariantCaller

# build() only checks that input files exist, so empty placeholder files are enough.


@pytest.fixture
def files(tmp_path):
    names = ["ref.fasta", "sample.bam", "sample.bam.bai", "sample_R1.fastq.gz", "sample_R2.fastq.gz"]
    for name in names:
        (tmp_path / name).touch()
    return {name: str(tmp_path / name) for name in names}


def fastqs(files):
    return [
        FastqReader(files["sample_R1.fastq.gz"], read=1),
        FastqReader(files["sample_R2.fastq.gz"], read=2),
    ]


def mapper(files, read_groups={"SM": "sample"}):
    return Mapper(library="bwa", input_step=fastqs(files), params={"read_groups": read_groups})


@pytest.mark.parametrize("library", ["deepvariant", "haplotypecaller"])
def test_bam_pipeline_builds(files, tmp_path, library):
    caller = VariantCaller(library=library, normal_sample=BamReader(files["sample.bam"]))
    params = Pipeline(ref_fasta=files["ref.fasta"]).add(caller).build(workdir=str(tmp_path / "work"))
    assert (tmp_path / "work" / "params.json").samefile(params)


def test_fastq_pipeline_builds(files, tmp_path):
    caller = VariantCaller(library="deepvariant", normal_sample=fastqs(files))
    Pipeline(ref_fasta=files["ref.fasta"]).add(caller).build(workdir=str(tmp_path / "work"))


def test_mapper_pipeline_builds(files, tmp_path):
    caller = VariantCaller(library="deepvariant", normal_sample=mapper(files))
    Pipeline(ref_fasta=files["ref.fasta"]).add(caller).build(workdir=str(tmp_path / "work"))


def test_alignment_only_builds(files, tmp_path):
    Pipeline(ref_fasta=files["ref.fasta"]).add(mapper(files)).build(workdir=str(tmp_path / "work"))


def test_unknown_caller(files):
    with pytest.raises(ValueError):
        VariantCaller(library="nosuchcaller", normal_sample=BamReader(files["sample.bam"]))


def test_unknown_mapper(files):
    with pytest.raises(ValueError):
        Mapper(library="nosuchmapper", input_step=fastqs(files), params={"read_groups": {"SM": "s"}})


def test_mapper_needs_sample_name(files):
    with pytest.raises(ValueError):
        mapper(files, read_groups={"ID": "run1"})


def test_fastq_reads_must_be_1_and_2(files):
    reads = [
        FastqReader(files["sample_R1.fastq.gz"], read=1),
        FastqReader(files["sample_R2.fastq.gz"], read=1),
    ]
    with pytest.raises(ValueError):
        VariantCaller(library="deepvariant", normal_sample=reads)


def test_haplotypecaller_rejects_model_type(files):
    with pytest.raises(ValueError):
        VariantCaller(
            library="haplotypecaller",
            normal_sample=BamReader(files["sample.bam"]),
            params={"model_type": "WES"},
        )


def test_missing_bam_index(files, tmp_path):
    (tmp_path / "sample.bam.bai").unlink()
    caller = VariantCaller(library="deepvariant", normal_sample=BamReader(files["sample.bam"]))
    with pytest.raises(FileNotFoundError):
        Pipeline(ref_fasta=files["ref.fasta"]).add(caller).build(workdir=str(tmp_path / "work"))


def test_missing_reference(files, tmp_path):
    caller = VariantCaller(library="deepvariant", normal_sample=BamReader(files["sample.bam"]))
    with pytest.raises(FileNotFoundError):
        Pipeline(ref_fasta=str(tmp_path / "nope.fasta")).add(caller).build(workdir=str(tmp_path / "work"))
