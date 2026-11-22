# COSAP-NX API Reference

Everything is imported from `cosap_nx`:

```python
from cosap_nx import (
    FastqReader, BamReader,
    Trimmer, Mapper, MDUP, Recalibrator, VariantCaller,
    Pipeline, PipelineRunner, VariantComparator,
)
```

Steps point at the step(s) they read from, like in COSAP. Nothing runs until `PipelineRunner.run_pipeline()`.

---

## Inputs

### FastqReader

```python
FastqReader(filename: str, read: int, name: str = None, platform: str = "ILLUMINA")
```

One FASTQ file. Paired-end input is a list of two readers with `read=1` and `read=2`.

| Argument | Description |
|----------|-------------|
| `filename` | Path to the FASTQ file. Must exist. |
| `read` | 1 or 2. |
| `name` | Sample name. If not given, taken from the filename: `sample_R1.fastq.gz` → `sample`. Give both reads the same name. |

### BamReader

```python
BamReader(filename: str, name: str = None)
```

An aligned BAM. `name` defaults to the filename without `.bam`. `build()` checks that the BAM and its index (`<bam>.bai` or `<name>.bai`) exist.

---

## Steps

Every step has a `name`. If you don't give one, it is built from the input's name like in COSAP. Names only have to be unique within a stage, but two different steps with the same name in one stage is an error.

### Trimmer

```python
Trimmer(input_step: List[FastqReader], name: str = None)
```

Read trimming (fastp in COSAP). Default name: the sample name. Not supported in real runs yet.

### Mapper

```python
Mapper(library: str, input_step: Union[List[FastqReader], Trimmer], params: dict = {}, name: str = None)
```

| Argument | Description |
|----------|-------------|
| `library` | `"bwa"` or `"bwa2"`. Other names (e.g. `"bowtie2"`) only run with `-stub-run`. |
| `input_step` | A FASTQ pair or a `Trimmer`. |
| `params` | Must contain `read_groups` with at least `SM`, like in COSAP. Other read group fields (`ID`, `LB`, `PL`, `PU`) are optional. |

Default name: `<input>_<library>`, e.g. `sample_bwa`.

### MDUP

```python
MDUP(input_step: Union[Mapper, BamReader], name: str = None, spark: bool = True)
```

Duplicate marking. `spark=True` is GATK MarkDuplicatesSpark, `False` is Picard MarkDuplicates, as in COSAP. Keeps its input's name. Not supported in real runs yet.

### Recalibrator

```python
Recalibrator(input_step: Union[MDUP, Mapper, BamReader], name: str = None)
```

Base quality score recalibration (GATK BQSR). Keeps its input's name. Not supported in real runs yet.

### VariantCaller

```python
VariantCaller(
    library: str,
    germline: Union[Mapper, MDUP, Recalibrator, BamReader] = None,
    tumor: Union[Mapper, MDUP, Recalibrator, BamReader] = None,
    params: dict = {},
    name: str = None,
    gvcf: bool = False,
)
```

| Argument | Description |
|----------|-------------|
| `library` | `"deepvariant"` or `"haplotypecaller"`. Others (e.g. `"mutect"`, `"strelka"`) only run with `-stub-run`. |
| `germline` | The normal sample. Germline callers only use this one. |
| `tumor` | The tumor sample, for somatic callers. |
| `params` | DeepVariant: `model_type` (`WGS`, `WES` or `PACBIO`, default `WGS`). |
| `gvcf` | DeepVariant writes a gVCF next to the VCF. HaplotypeCaller writes only the gVCF, as in COSAP. |

At least one of `germline` and `tumor` is required. Default name: the input names and the library joined, e.g. `sample_bwa_deepvariant` or `normal_bwa_tumor_bwa_mutect`.

---

## Pipeline

```python
Pipeline(ref_fasta: str = None)
```

### `add(step) -> Pipeline`

Adds a step. Steps it reads from are included automatically, so adding only the callers is enough. Returns the pipeline, so calls can be chained.

### `set_reference(ref_fasta) -> Pipeline`

Sets the reference fasta.

### `build(workdir: str, cpus: int = 4, memory: str = "8 GB") -> str`

Checks the inputs, `<ref>` and `<ref>.fai`, then writes `<workdir>/params.json` and returns its path. `cpus` and `memory` are used for every step. Nextflow's local executor refuses tasks that ask for more than the machine has.

```python
pipeline = Pipeline(ref_fasta="ref.fasta")
for m in ["bwa", "bwa2"]:
    mapped = Mapper(library=m, input_step=reads, params={"read_groups": {"SM": "sample"}})
    for c in ["deepvariant", "haplotypecaller"]:
        pipeline.add(VariantCaller(library=c, germline=mapped))
config = pipeline.build(workdir="./workdir")
```

---

## PipelineRunner

```python
PipelineRunner(nextflow_path: str = None)
```

Finds `nextflow` on the PATH if no path is given (raises `NextflowNotFoundError` otherwise).

### `run_pipeline(pipeline_config, profile="docker", resume=True, work_dir=None, compare=False, truth_vcf=None) -> Path`

Runs `nextflow run nf/main.nf -params-file <config>` in the folder of the config file and returns the `results/` directory. Raises `PipelineExecutionError` if Nextflow fails.

| Argument | Description |
|----------|-------------|
| `profile` | `docker` or `singularity`. |
| `resume` | Pass `-resume` to reuse finished tasks. |
| `work_dir` | Nextflow work directory (default `work/` next to the config). |
| `compare` | Afterwards, compare each sample's VCFs into `results/comparison/<sample>/`. |
| `truth_vcf` | With `compare`, also compute precision/recall/F1 against this VCF. |

### `compare_results(results_dir, truth_vcf=None) -> Path`

The comparison step on its own, for a run that already finished. Writes `variant_counts.tsv`, `overlap_summary.json` and (with a truth VCF) `metrics_vs_truth.json` per sample.

---

## VariantComparator

```python
VariantComparator(variant_files: List[dict], bed_file: str = None)
```

Same input as COSAP's comparator: a list of `{"path", "caller", "mapper"}`. Each VCF is labelled `<mapper>_<caller>`. `bed_file` is accepted but not used.

### `VariantComparator.from_results(results_dir, sample=None)`

Builds the list from `results/variants.json`. With `sample`, only that sample's VCFs are used. Raises `ValueError` if two VCFs would get the same label (e.g. two samples both have `bwa_deepvariant`); pick a sample then. VCFs that came from a BAM input get the mapper label `bam`.

### Metrics

| Method | Output |
|--------|--------|
| `compute_statistics(output_dir)` | `variant_counts.tsv`: SNPs, indels, totals |
| `compute_overlap(output_dir)` | `overlap_summary.json`: pairwise shared/unique counts and Jaccard |
| `compute_metrics_vs_truth(truth_vcf, output_dir)` | `metrics_vs_truth.json`: TP, FP, FN, precision, recall, F1. Handles `chr20` vs `20` naming. |

### Plots

`draw_venn_diagram`, `draw_venn2_plot`, `draw_venn3_plot`, `draw_upset_plot`, `draw_similarity_plot`, `draw_precision_recall_plot(truth_vcf)`. Each takes an `output_file` (PNG).

### Other

- `get_variants()`: all variants as a DataFrame.
- `create_intersection_bed(query, output_file)`: export a set of variants as BED, e.g. `"bwa_deepvariant & ~bwa_haplotypecaller"`.

---

## params.json

What `build()` writes and `nf/main.nf` reads:

```json
{
  "version": "0.3.0",
  "outdir": "/path/workdir/results",
  "defaults": {"cpus": 4, "memory": "8 GB"},
  "references": {"fasta": "/refs/ref.fasta"},
  "inputs": {
    "sample": {"type": "fastq", "r1": "...", "r2": "..."}
  },
  "mapping": {
    "sample_bwa": {"input": "inputs/sample", "library": "bwa", "params": {"read_groups": {"SM": "sample"}}}
  },
  "variant_calling": {
    "sample_bwa_deepvariant": {"normal": "mapping/sample_bwa", "library": "deepvariant", "params": {}}
  }
}
```

Each stage (`trim`, `mapping`, `mdup`, `calibrate`, `variant_calling`) maps step names to their settings. Steps refer to their inputs as `"stage/name"`. A step can also have its own `cpus` and `memory`, which override `defaults`.

---

## CLI

```bash
cosap-nx-run --config params.json [--profile docker] [--no-resume] [--work-dir DIR] [--compare] [--truth-vcf VCF]
```
