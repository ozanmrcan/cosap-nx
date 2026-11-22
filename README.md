# COSAP-NX v0.3.0

A COSAP-like Python API for building NGS pipelines, with Nextflow running them instead of Snakemake.

You describe the pipeline in Python the same way as in [COSAP](https://github.com/MBaysanLab/cosap): every step points at the step(s) it reads from. `Pipeline.build()` writes a `params.json`, and `PipelineRunner` runs `nf/main.nf` on it with each tool in its own container.

## Requirements

- Python 3.8+
- [Nextflow](https://www.nextflow.io/) 23.04 or newer (Java 11+)
- Docker

## Installation

```bash
cd cosap-nx
python3 -m venv venv
source venv/bin/activate
pip install -e .
```

See [docs/SETUP.md](docs/SETUP.md) for more detail.

## Quick start

Two mappers × two callers on one FASTQ pair, compared at the end:

```python
from cosap_nx import FastqReader, Mapper, VariantCaller, Pipeline, PipelineRunner

reads = [
    FastqReader("/data/sample_R1.fastq.gz", read=1, name="sample"),
    FastqReader("/data/sample_R2.fastq.gz", read=2, name="sample"),
]

pipeline = Pipeline(ref_fasta="/refs/ref.fasta")
for mapper in ["bwa", "bwa2"]:
    mapped = Mapper(library=mapper, input_step=reads, params={"read_groups": {"SM": "sample"}})
    for caller in ["deepvariant", "haplotypecaller"]:
        pipeline.add(VariantCaller(library=caller, germline=mapped))

config = pipeline.build(workdir="./workdir")
results = PipelineRunner().run_pipeline(config, compare=True, truth_vcf="/data/truth.vcf.gz")
```

The mapper is added automatically because the callers read from it. Steps get COSAP-style names (`sample_bwa`, `sample_bwa_deepvariant`), and each mapper runs once even though two callers read from it.

Starting from a BAM:

```python
from cosap_nx import BamReader, VariantCaller, Pipeline, PipelineRunner

bam = BamReader("/data/NA12878.bam", name="NA12878")
caller = VariantCaller(library="deepvariant", germline=bam, params={"model_type": "WGS"})

config = Pipeline(ref_fasta="/refs/ref.fasta").add(caller).build(workdir="./workdir")
PipelineRunner().run_pipeline(config)
```

More in [examples/](examples/) and the [user guide](docs/USER_GUIDE.md).

## Steps and tools

| Step | Stage | Tools that run for real |
|------|-------|-------------------------|
| `Trimmer(input_step)` | trim | (none yet) |
| `Mapper(library, input_step, params)` | mapping | `bwa`, `bwa2` |
| `MDUP(input_step, spark=True)` | mdup | (none yet) |
| `Recalibrator(input_step)` | calibrate | (none yet) |
| `VariantCaller(library, germline, tumor, params, gvcf)` | variant_calling | `deepvariant`, `haplotypecaller` |

Any other step or tool can already be put in a pipeline and tested with `-stub-run`, which runs placeholders instead of tools. A normal run stops and tells you which tool isn't supported yet.

## Comparing results

After variant calling, the workflow writes `results/variants.json`, which lists every VCF with its mapper, caller and sample. With `compare=True` (or `cosap-nx-run --compare`) the runner compares each sample's VCFs into `results/comparison/<sample>/`:

- `variant_counts.tsv`: SNPs, indels and totals per pipeline
- `overlap_summary.json`: pairwise overlap and Jaccard similarity
- `metrics_vs_truth.json`: precision, recall and F1, if a truth VCF is given

For plots or anything else, build the comparator yourself:

```python
from cosap_nx import VariantComparator

comparator = VariantComparator.from_results(results, sample="sample")
comparator.draw_upset_plot(output_file="upset.png")
comparator.draw_venn_diagram(output_file="venn.png")
```

`VariantComparator(variant_files=[{"path", "caller", "mapper"}, ...])` also works the way it does in COSAP.

## Output

```
workdir/
  params.json
  results/
    bam/<mapper>/<step>.bam(.bai)
    vcf/<caller>/<step>.vcf.gz(.tbi)      # .g.vcf.gz with gvcf=True
    variants.json
    comparison/<sample>/                  # with compare=True
    logs/trace.txt, report.html, timeline.html
```

## Other ways to run

```bash
# CLI runner on an existing params.json
cosap-nx-run --config workdir/params.json --profile docker --compare --truth-vcf truth.vcf.gz

# Nextflow directly
nextflow run nf/main.nf -params-file examples/sample_params.json -profile docker

# Wiring only: no tools or containers
nextflow run nf/main.nf -params-file workdir/params.json -stub-run
```

## Reference files

- `<ref>.fai` next to the fasta.
- `bwa`: BWA index files (`bwa index <ref>`).
- `bwa2`: BWA-MEM2 index files (`bwa-mem2 index <ref>`).
- `haplotypecaller`: a sequence dictionary (`<ref>.dict` or `<base>.dict`).
- BAM inputs need a `.bai`.

## Tests

```bash
pip install -e ".[dev]"
pytest              # API and comparator, about a second
pytest -m stub      # workflow wiring with -stub-run, about 20 s
pytest -m slow      # one real chr20 run, about 10 min, needs Docker and ../data
```

## Project structure

```
cosap-nx/
  cosap_nx/
    api.py           # step classes and Pipeline
    runner.py        # PipelineRunner, cosap-nx-run
    comparator.py    # VariantComparator
  nf/
    main.nf          # the workflow
    modules/         # one process per tool
    nextflow.config
    conf/docker.config
  examples/
  tests/
```

## Documentation

- [Setup](docs/SETUP.md)
- [User guide](docs/USER_GUIDE.md)
- [API reference](docs/API_REFERENCE.md)

## Version history

- **v0.3.0**: COSAP-style step graph (Trimmer, Mapper, MDUP, Recalibrator, VariantCaller with germline/tumor), many steps per stage in one run, BWA-MEM2, automatic comparison after a run.
- **v0.2.0**: VariantComparator and HaplotypeCaller.
- **v0.1.2**: FASTQ input with BWA-MEM.
- **v0.1.0**: BAM → DeepVariant → VCF.

See [CHANGELOG.md](CHANGELOG.md).
