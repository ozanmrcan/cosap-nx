# COSAP-NX User Guide

This guide explains how to use COSAP-NX to run variant calling pipelines.

---

## Overview

COSAP-NX provides a **COSAP-like Python API** that generates and runs **Nextflow** pipelines under the hood.

```
Your Python Script          COSAP-NX              Nextflow           Docker
┌─────────────────┐    ┌──────────────┐    ┌──────────────┐    ┌─────────────┐
│ BamReader       │───▶│ Generates    │───▶│ Executes     │───▶│ DeepVariant │
│ VariantCaller   │    │ params.json  │    │ main.nf      │    │ Container   │
│ Pipeline        │    │              │    │              │    │             │
└─────────────────┘    └──────────────┘    └──────────────┘    └─────────────┘
```

---

## Quick Start

### 1. Activate the Virtual Environment

Every time you open a new terminal:

```bash
cd ~/grad_project/cosap-nx
source venv/bin/activate
```

You'll see `(venv)` in your prompt when active.

### 2. Basic Pipeline Script

```python
from cosap_nx import BamReader, VariantCaller, Pipeline, PipelineRunner

# 1. Define your input BAM file
bam = BamReader("/path/to/sample.bam")

# 2. Configure the variant caller
caller = VariantCaller(
    library="deepvariant",
    normal_sample=bam
)

# 3. Build the pipeline
pipeline = Pipeline()
pipeline.set_reference("/path/to/reference.fa")
pipeline.add(caller)

# 4. Generate config and run
config_path = pipeline.build(workdir="./my_analysis")

runner = PipelineRunner()
runner.run_pipeline(config_path, profile="docker")
```

---

## Step-by-Step Explanation

### Step 1: BamReader

`BamReader` wraps an input BAM file:

```python
from cosap_nx import BamReader

# Basic usage - name is auto-extracted from filename
bam = BamReader("/data/NA12878.bam")
print(bam.name)  # "NA12878"

# Custom name
bam = BamReader("/data/sample.bam", name="Patient_001")
print(bam.name)  # "Patient_001"
```

**Requirements:**
- BAM file must exist
- BAM index must exist (`.bam.bai` or `.bai`)

### Step 1b: FastqReader (v0.1.2+)

**NEW in v0.1.2:** You can now start from FASTQ files instead of BAM!

`FastqReader` wraps FASTQ files for alignment:

```python
from cosap_nx import FastqReader

# Paired-end reads (most common)
fastq_r1 = FastqReader("/data/sample_R1.fastq.gz", read=1, name="MySample")
fastq_r2 = FastqReader("/data/sample_R2.fastq.gz", read=2, name="MySample")

# Name is auto-extracted if not provided
fastq_r1 = FastqReader("/data/NA12878_R1.fastq.gz", read=1)
print(fastq_r1.name)  # "NA12878"
```

**Key points:**
- Always create **2 FastqReader instances** for paired-end data (read=1 and read=2)
- Both must have the **same name** (will be checked during validation)
- FASTQ files can be gzipped (`.fastq.gz`) or uncompressed (`.fastq`)

**Prerequisites for FASTQ mode:**
The reference genome must be indexed with BWA:

```bash
bwa index /path/to/reference.fasta
```

This creates `.amb`, `.ann`, `.bwt`, `.pac`, `.sa` index files.

**Example: FASTQ to VCF**

```python
from cosap_nx import FastqReader, VariantCaller, Pipeline, PipelineRunner

# Define paired-end FASTQ inputs
fastq_r1 = FastqReader("sample_R1.fastq.gz", read=1, name="MySample")
fastq_r2 = FastqReader("sample_R2.fastq.gz", read=2, name="MySample")

# Pass list of FastqReaders to VariantCaller
caller = VariantCaller(
    library="deepvariant",
    normal_sample=[fastq_r1, fastq_r2],  # List, not single file!
    params={"model_type": "WGS"}
)

# Build and run (same as BAM mode)
pipeline = Pipeline(ref_fasta="reference.fasta")
pipeline.add(caller)
config_path = pipeline.build(workdir="./workdir", cpus=8, memory="16 GB")

runner = PipelineRunner()
runner.run_pipeline(config_path)
```

**Outputs:**
- Aligned BAM: `workdir/results/bam/MySample.bam`
- Variants VCF: `workdir/results/vcf/deepvariant/MySample.vcf.gz`

### Step 2: VariantCaller

`VariantCaller` configures which tool to use and how.

**Accepts either:**
- **BAM input:** Single `BamReader` instance (v0.1.0 behavior)
- **FASTQ input:** List of 2 `FastqReader` instances (v0.1.2+)

```python
from cosap_nx import VariantCaller

# Option 1: BAM input
caller = VariantCaller(
    library="deepvariant",
    normal_sample=bam  # Single BamReader
)

# Option 2: FASTQ input (v0.1.2+)
caller = VariantCaller(
    library="deepvariant",
    normal_sample=[fastq_r1, fastq_r2]  # List of FastqReaders
)

# With custom parameters
caller = VariantCaller(
    library="deepvariant",
    normal_sample=bam,  # or [fastq_r1, fastq_r2]
    params={
        "model_type": "WES",           # WGS, WES, or PACBIO
        "germline_sample_name": "NA12878"
    },
    gvcf=True  # Also output gVCF
)
```

**Supported libraries (v0.1.2):**
- `deepvariant` - Google's deep learning variant caller

**Model types:**
| Type | Use for |
|------|---------|
| `WGS` | Whole genome sequencing (default) |
| `WES` | Whole exome sequencing |
| `PACBIO` | PacBio long reads |

### Step 3: Pipeline

`Pipeline` collects steps and builds the configuration:

```python
from cosap_nx import Pipeline

pipeline = Pipeline()

# Set reference genome (required)
pipeline.set_reference("/refs/hg38.fa")

# Add variant caller step
pipeline.add(caller)

# Build configuration
config_path = pipeline.build(
    workdir="./my_analysis",  # Where outputs go
    cpus=8,                   # CPUs for DeepVariant
    memory="16 GB"            # Memory allocation
)
```

**Reference requirements:**
- FASTA file (`.fa` or `.fasta`)
- Index file (`.fai`) - create with `samtools faidx ref.fa`

### Step 4: PipelineRunner

`PipelineRunner` executes the Nextflow workflow:

```python
from cosap_nx import PipelineRunner

runner = PipelineRunner()

# Run with Docker (default)
output_dir = runner.run_pipeline(
    pipeline_config=config_path,
    profile="docker"
)

# Other options
output_dir = runner.run_pipeline(
    pipeline_config=config_path,
    profile="docker",
    resume=True,       # Resume from cache (default: True)
    work_dir="./work"  # Nextflow work directory
)
```

**Profiles:**
- `docker` - Run tools in Docker containers (recommended)
- `singularity` - Run tools in Singularity containers (for HPC)

---

## Output Structure

After a successful run:

```
my_analysis/
├── params.json                    # Generated configuration
├── results/
│   ├── vcf/
│   │   └── deepvariant/
│   │       ├── NA12878.vcf.gz     # Variant calls
│   │       └── NA12878.vcf.gz.tbi # VCF index
│   └── logs/
│       ├── trace.txt              # Execution trace
│       ├── report.html            # Visual report
│       └── timeline.html          # Timeline visualization
└── work/                          # Nextflow work directory (can be deleted)
```

---

## Complete Examples

### Example 1: Whole Genome Sequencing

```python
from cosap_nx import BamReader, VariantCaller, Pipeline, PipelineRunner

bam = BamReader("/data/WGS_sample.bam")

caller = VariantCaller(
    library="deepvariant",
    normal_sample=bam,
    params={"model_type": "WGS"}
)

pipeline = Pipeline()
pipeline.set_reference("/refs/GRCh38.fa")
pipeline.add(caller)

config = pipeline.build(
    workdir="./wgs_analysis",
    cpus=16,
    memory="32 GB"
)

runner = PipelineRunner()
runner.run_pipeline(config, profile="docker")
```

### Example 2: Whole Exome Sequencing

```python
from cosap_nx import BamReader, VariantCaller, Pipeline, PipelineRunner

bam = BamReader("/data/exome_sample.bam", name="Patient_42")

caller = VariantCaller(
    library="deepvariant",
    normal_sample=bam,
    params={"model_type": "WES"},
    gvcf=True  # Also produce gVCF for joint calling
)

pipeline = Pipeline()
pipeline.set_reference("/refs/hg38.fa")
pipeline.add(caller)

config = pipeline.build(workdir="./exome_analysis")

runner = PipelineRunner()
runner.run_pipeline(config)
```

### Example 3: Using the CLI Runner

Instead of writing Python, you can run directly with a JSON config:

```bash
# Create params.json manually or from a previous run
cosap-nx-run --config params.json --profile docker
```

---

## Running Nextflow Directly (Advanced)

For debugging or manual testing, you can run Nextflow without the Python wrapper:

```bash
cd ~/grad_project/cosap-nx

nextflow run nf/main.nf \
    --bam /path/to/sample.bam \
    --ref_fasta /path/to/reference.fa \
    --sample_id my_sample \
    --outdir ./results \
    --model_type WGS \
    -profile docker
```

Or with a params file:

```bash
nextflow run nf/main.nf \
    -params-file params.json \
    -profile docker
```

---

## Tips and Best Practices

### 1. Resource Allocation

DeepVariant is resource-intensive. Recommended settings:

| Data Type | CPUs | Memory |
|-----------|------|--------|
| WES | 8 | 16 GB |
| WGS (30x) | 16+ | 32+ GB |

### 2. Resuming Failed Runs

Nextflow caches intermediate results. If a run fails, just run again - it will resume from where it left off:

```python
runner.run_pipeline(config, resume=True)  # default
```

### 3. Cleaning Up

The `work/` directory can get large. After a successful run:

```bash
rm -rf work/
```

### 4. Checking Logs

If something fails, check:
1. Terminal output for error messages
2. `results/logs/trace.txt` for detailed execution trace
3. `.nextflow.log` for Nextflow-specific logs

---

## Next Steps

- [API Reference](API_REFERENCE.md) - Detailed class and method documentation
- [SETUP.md](SETUP.md) - Installation instructions
