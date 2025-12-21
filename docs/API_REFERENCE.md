# COSAP-NX API Reference

Complete reference for all COSAP-NX classes and methods.

---

## Module: cosap_nx

```python
from cosap_nx import BamReader, VariantCaller, Pipeline, PipelineRunner
```

---

## BamReader

Represents an input BAM file.

### Constructor

```python
BamReader(filename: str, name: str = None)
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `filename` | `str` | Yes | Path to BAM file |
| `name` | `str` | No | Sample name. If not provided, extracted from filename |

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `filename` | `str` | Absolute path to BAM file |
| `name` | `str` | Sample identifier |

### Methods

#### `get_output() -> str`

Returns the absolute path to the BAM file.

```python
bam = BamReader("/data/sample.bam")
print(bam.get_output())  # "/data/sample.bam"
```

#### `validate() -> None`

Checks that the BAM file and its index exist. Raises `FileNotFoundError` if not found.

```python
bam = BamReader("/data/sample.bam")
bam.validate()  # Raises FileNotFoundError if file or index missing
```

### Example

```python
from cosap_nx import BamReader

# Auto-extract name from filename
bam1 = BamReader("/data/NA12878_sorted.bam")
print(bam1.name)  # "NA12878_sorted"

# Specify custom name
bam2 = BamReader("/data/sample.bam", name="Patient_001")
print(bam2.name)  # "Patient_001"
```

---

## VariantCaller

Configures variant calling for a sample.

### Constructor

```python
VariantCaller(
    library: str,
    normal_sample: BamReader = None,
    germline: BamReader = None,      # Alias for normal_sample
    tumor_sample: BamReader = None,  # Not supported in v0.1
    params: Dict[str, Any] = {},
    gvcf: bool = False,
    name: str = None
)
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `library` | `str` | Yes | Variant caller to use. v0.1 supports: `"deepvariant"` |
| `normal_sample` | `BamReader` | Yes | Input BAM for germline calling |
| `germline` | `BamReader` | No | Alias for `normal_sample` (COSAP compatibility) |
| `tumor_sample` | `BamReader` | No | Not supported in v0.1 |
| `params` | `dict` | No | Tool-specific parameters |
| `gvcf` | `bool` | No | Output gVCF in addition to VCF (default: `False`) |
| `name` | `str` | No | Step name. Auto-generated if not provided |

### Supported Parameters

For `library="deepvariant"`:

| Parameter | Type | Values | Default | Description |
|-----------|------|--------|---------|-------------|
| `model_type` | `str` | `"WGS"`, `"WES"`, `"PACBIO"` | `"WGS"` | Sequencing type |
| `germline_sample_name` | `str` | any | from BAM | Sample name in output VCF |

### Methods

#### `get_sample_id() -> str`

Returns the sample identifier for the output VCF.

#### `get_model_type() -> str`

Returns the DeepVariant model type.

### Example

```python
from cosap_nx import BamReader, VariantCaller

bam = BamReader("/data/exome.bam")

# Basic usage
caller = VariantCaller(
    library="deepvariant",
    normal_sample=bam
)

# With all options
caller = VariantCaller(
    library="deepvariant",
    normal_sample=bam,
    params={
        "model_type": "WES",
        "germline_sample_name": "Sample_001"
    },
    gvcf=True,
    name="exome_deepvariant"
)
```

---

## Pipeline

Collects pipeline steps and builds configuration.

### Constructor

```python
Pipeline(ref_fasta: str = None)
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `ref_fasta` | `str` | No | Path to reference FASTA (can also use `set_reference()`) |

### Methods

#### `set_reference(ref_fasta: str) -> Pipeline`

Set the reference genome path. Returns self for method chaining.

```python
pipeline = Pipeline()
pipeline.set_reference("/refs/hg38.fa")

# Or chain
pipeline = Pipeline().set_reference("/refs/hg38.fa")
```

#### `add(step: VariantCaller) -> Pipeline`

Add a step to the pipeline. Returns self for method chaining.

```python
pipeline.add(caller)

# Or chain
pipeline = Pipeline().set_reference("/refs/hg38.fa").add(caller)
```

#### `build(workdir: str, cpus: int = 8, memory: str = "16 GB") -> str`

Generate the pipeline configuration file. Returns path to `params.json`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `workdir` | `str` | required | Output directory |
| `cpus` | `int` | `8` | CPU cores for variant calling |
| `memory` | `str` | `"16 GB"` | Memory allocation |

```python
config_path = pipeline.build(
    workdir="./analysis",
    cpus=16,
    memory="32 GB"
)
# Returns: "./analysis/params.json"
```

#### `get_workdir() -> Optional[str]`

Returns the workdir path set during `build()`, or `None` if not yet built.

### Example

```python
from cosap_nx import BamReader, VariantCaller, Pipeline

bam = BamReader("/data/sample.bam")
caller = VariantCaller(library="deepvariant", normal_sample=bam)

# Method 1: Step by step
pipeline = Pipeline()
pipeline.set_reference("/refs/hg38.fa")
pipeline.add(caller)
config = pipeline.build(workdir="./output")

# Method 2: Chained
config = (
    Pipeline()
    .set_reference("/refs/hg38.fa")
    .add(caller)
    .build(workdir="./output")
)
```

---

## PipelineRunner

Executes Nextflow pipelines.

### Constructor

```python
PipelineRunner(nextflow_path: str = None)
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `nextflow_path` | `str` | No | Path to Nextflow executable. Auto-detected if not provided |

### Methods

#### `run_pipeline(pipeline_config: str, profile: str = "docker", resume: bool = True, work_dir: str = None) -> Path`

Execute the pipeline. Returns path to output directory.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pipeline_config` | `str` | required | Path to `params.json` |
| `profile` | `str` | `"docker"` | Nextflow profile (`"docker"`, `"singularity"`) |
| `resume` | `bool` | `True` | Resume from cached results |
| `work_dir` | `str` | `None` | Nextflow work directory |

### Exceptions

| Exception | When |
|-----------|------|
| `NextflowNotFoundError` | Nextflow not installed or not in PATH |
| `FileNotFoundError` | Config file doesn't exist |
| `PipelineExecutionError` | Pipeline fails during execution |

### Example

```python
from cosap_nx import PipelineRunner

runner = PipelineRunner()

# Basic run
output_dir = runner.run_pipeline("./analysis/params.json")

# With options
output_dir = runner.run_pipeline(
    pipeline_config="./analysis/params.json",
    profile="docker",
    resume=True,
    work_dir="./custom_work"
)

print(f"VCF files in: {output_dir}")
```

---

## Exceptions

### NextflowNotFoundError

Raised when Nextflow executable cannot be found.

```python
from cosap_nx.runner import NextflowNotFoundError

try:
    runner = PipelineRunner()
except NextflowNotFoundError:
    print("Please install Nextflow")
```

### PipelineExecutionError

Raised when the Nextflow pipeline fails.

```python
from cosap_nx.runner import PipelineExecutionError

try:
    runner.run_pipeline(config)
except PipelineExecutionError as e:
    print(f"Pipeline failed: {e}")
```

---

## Generated Configuration (params.json)

The `Pipeline.build()` method generates a JSON configuration:

```json
{
  "sample_id": "NA12878",
  "bam": "/data/NA12878.bam",
  "ref_fasta": "/refs/hg38.fa",
  "outdir": "/analysis/results",
  "model_type": "WGS",
  "gvcf": false,
  "cpus": 8,
  "memory": "16 GB"
}
```

This file can also be created manually and used with the CLI or direct Nextflow execution.

---

## CLI Interface

### cosap-nx-run

Command-line tool for running pipelines.

```bash
cosap-nx-run --config <path> [options]
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--config` | `-c` | required | Path to params.json |
| `--profile` | `-p` | `docker` | Nextflow profile |
| `--no-resume` | | | Disable resume from cache |
| `--work-dir` | `-w` | | Custom work directory |

```bash
# Examples
cosap-nx-run -c params.json
cosap-nx-run -c params.json -p singularity
cosap-nx-run -c params.json --no-resume -w /tmp/work
```
