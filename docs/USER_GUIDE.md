# COSAP-NX User Guide

How to build and run pipelines with COSAP-NX. For every argument, see the [API reference](API_REFERENCE.md).

---

## How it works

```
your script  →  params.json  →  nextflow run nf/main.nf  →  tools in containers  →  workdir/results/
```

1. In Python you create steps. Each step points at the step(s) it reads from, like in COSAP.
2. `Pipeline.build()` follows those links back to the inputs and writes `params.json`, which lists every step by stage.
3. `PipelineRunner.run_pipeline()` runs Nextflow on that file. Stages run in a fixed order (trim → mapping → mdup → calibrate → variant_calling), and each stage can hold any number of steps.

Python doesn't know which tools exist. The workflow checks that before running anything and stops with a clear error if a tool isn't supported yet.

Before you start, activate the venv in `cosap-nx/`:

```bash
source venv/bin/activate
```

---

## Inputs

FASTQ pairs are a list of two `FastqReader`s:

```python
from cosap_nx import FastqReader, BamReader

reads = [
    FastqReader("/data/sample_R1.fastq.gz", read=1, name="sample"),
    FastqReader("/data/sample_R2.fastq.gz", read=2, name="sample"),
]
```

A BAM needs its index next to it:

```python
bam = BamReader("/data/sample.bam", name="sample")
```

The names are used to name everything downstream, so keep them short and unique. If a BAM and a FASTQ pair are the same sample, give them different names (e.g. `sample` and `sample_bam`). Otherwise `build()` stops because two inputs share a name.

---

## A single caller

```python
from cosap_nx import BamReader, VariantCaller, Pipeline, PipelineRunner

bam = BamReader("/data/sample.bam", name="sample")
caller = VariantCaller(library="deepvariant", germline=bam, params={"model_type": "WES"})

config = Pipeline(ref_fasta="/refs/ref.fasta").add(caller).build(workdir="./analysis")
PipelineRunner().run_pipeline(config)
```

The VCF ends up in `analysis/results/vcf/deepvariant/sample_deepvariant.vcf.gz`.

---

## From FASTQ

A `Mapper` needs a sample name in its read groups, as in COSAP:

```python
from cosap_nx import Mapper, VariantCaller, Pipeline

mapped = Mapper(library="bwa", input_step=reads, params={"read_groups": {"SM": "sample"}})
caller = VariantCaller(library="haplotypecaller", germline=mapped)

config = Pipeline(ref_fasta="/refs/ref.fasta").add(caller).build(workdir="./analysis")
```

You only need to add the caller. The mapper comes with it because the caller reads from it. The BAM is kept in `results/bam/bwa/sample_bwa.bam`.

The reference needs a BWA index for `bwa` (`bwa index ref.fasta`) and a BWA-MEM2 index for `bwa2` (`bwa-mem2 index ref.fasta`).

---

## Several tools in one run

This is the main difference from a fixed pipeline: any number of mappers and callers can run side by side. Expand the combinations in Python, the same way COSAP's `DNAPipeline` does:

```python
pipeline = Pipeline(ref_fasta="/refs/ref.fasta")
for m in ["bwa", "bwa2"]:
    mapped = Mapper(library=m, input_step=reads, params={"read_groups": {"SM": "sample"}})
    for c in ["deepvariant", "haplotypecaller"]:
        pipeline.add(VariantCaller(library=c, germline=mapped))

config = pipeline.build(workdir="./matrix")
```

This is one Nextflow run with four calling steps. Each mapper runs once, and its BAM goes to both callers.

---

## Tumor/normal

Somatic callers take both samples. The full COSAP preprocessing chain looks like this:

```python
from cosap_nx import Trimmer, Mapper, MDUP, Recalibrator, VariantCaller

def preprocess(reads, sample):
    mapped = Mapper(library="bwa", input_step=Trimmer(input_step=reads),
                    params={"read_groups": {"SM": sample}})
    return Recalibrator(input_step=MDUP(input_step=mapped))

mutect = VariantCaller(library="mutect", germline=preprocess(normal_reads, "normal"),
                       tumor=preprocess(tumor_reads, "tumor"))
```

The step waits until both samples are ready. Trimming, duplicate marking, BQSR and Mutect2 are not ported yet, so this pipeline only runs with `-stub-run` for now (see below).

---

## Comparing results

Pass `compare=True` to compare each sample's VCFs after the run:

```python
results = PipelineRunner().run_pipeline(config, compare=True, truth_vcf="/data/truth.vcf.gz")
```

This writes `results/comparison/<sample>/` with variant counts, pairwise overlap and Jaccard, and, with a truth VCF, precision/recall/F1. Pipelines are labelled `<mapper>_<caller>`, e.g. `bwa2_haplotypecaller`. Calls made straight from a BAM get the mapper label `bam`.

For plots, or to compare an older run, build the comparator from the results folder:

```python
from cosap_nx import VariantComparator

comparator = VariantComparator.from_results("./matrix/results", sample="sample")
comparator.draw_upset_plot(output_file="upset.png")
comparator.draw_similarity_plot(output_file="similarity.png")
comparator.draw_precision_recall_plot(truth_vcf="/data/truth.vcf.gz", output_file="pr.png")
```

You can also pass the VCFs by hand, as in COSAP: `VariantComparator(variant_files=[{"path": ..., "caller": ..., "mapper": ...}])`.

---

## Testing a pipeline without tools

`-stub-run` runs every step as a placeholder that only writes which inputs it got. This takes seconds and needs no containers, and it works for tools that aren't ported yet:

```bash
nextflow run nf/main.nf -params-file analysis/params.json -stub-run
```

Placeholder outputs go to `results/<stage>/`.

---

## Output

```
analysis/
  params.json
  work/                      Nextflow's task folders, can be deleted after a run
  results/
    bam/<mapper>/
    vcf/<caller>/
    variants.json            every VCF with mapper, caller and sample
    comparison/<sample>/     with compare=True
    logs/
      trace.txt              one line per task, with time, CPU and memory
      report.html
      timeline.html
```

---

## Tips

- **Resources.** `build(cpus=..., memory=...)` sets the default for every step. Nextflow's local executor won't start a task that asks for more than the machine has. To give one step more, add `cpus`/`memory` to that step in `params.json`.
- **Resuming.** `run_pipeline()` passes `-resume` by default, so after a failure only the unfinished steps run again.
- **When something fails,** check the terminal output, then `results/logs/trace.txt`, then `.nextflow.log` in the workdir.

---

## More

- [API reference](API_REFERENCE.md)
- [Setup](SETUP.md)
- `examples/` for runnable scripts on the chr20 test data
