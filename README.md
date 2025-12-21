# COSAP-NX v0.1

COSAP-like pipeline builder with **Nextflow** as the execution backend.

## v0.1 Scope

**BAM → DeepVariant → VCF** (germline variant calling)

This is a minimal vertical slice to prove the architecture before expanding.

## Requirements

- Python 3.8+
- [Nextflow](https://www.nextflow.io/) (>= 23.04)
- Docker (for `-profile docker`)
- Java 11+ (required by Nextflow)

## Installation

```bash
cd cosap-nx

# Create virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate

# Install COSAP-NX
pip install -e .
```

> **Note:** Always activate the venv (`source venv/bin/activate`) before working with COSAP-NX.
>
> See [docs/SETUP.md](docs/SETUP.md) for detailed installation instructions.

## Quick Start

### Option 1: Python API (COSAP-like)

```python
from cosap_nx import BamReader, VariantCaller, Pipeline, PipelineRunner

# Define input
bam = BamReader("/data/NA12878.bam", name="NA12878")

# Configure variant caller
dv = VariantCaller(
    library="deepvariant",
    normal_sample=bam,
    params={"model_type": "WGS"},
    gvcf=False
)

# Build and run pipeline
pipeline = Pipeline()
pipeline.set_reference("/refs/hg38.fa")
pipeline.add(dv)

config_path = pipeline.build(workdir="./workdir")

runner = PipelineRunner()
runner.run_pipeline(config_path, profile="docker")
```

### Option 2: Direct Nextflow (for testing)

```bash
nextflow run nf/main.nf \
    -params-file examples/sample_params.json \
    -profile docker
```

### Option 3: CLI Runner

```bash
cosap-nx-run --config params.json --profile docker
```

## Output Structure

```
workdir/
  params.json           # Generated config
  results/
    vcf/
      deepvariant/
        sample.vcf.gz
        sample.vcf.gz.tbi
    logs/
      trace.txt
      report.html
      timeline.html
```

## Configuration

### DeepVariant Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `model_type` | WGS, WES, or PACBIO | WGS |
| `gvcf` | Output gVCF format | false |
| `cpus` | CPU cores | 8 |
| `memory` | Memory allocation | 16 GB |

## Project Structure

```
cosap-nx/
  cosap_nx/
    __init__.py      # Public API exports
    api.py           # BamReader, VariantCaller, Pipeline
    runner.py        # PipelineRunner (calls Nextflow)
  nf/
    main.nf          # Nextflow workflow
    nextflow.config  # Root config
    conf/
      base.config    # Core settings
      docker.config  # Docker profile
  examples/
    germline_deepvariant.py
    sample_params.json
```

## Documentation

- [Setup Guide](docs/SETUP.md) - Installation and environment setup
- [User Guide](docs/USER_GUIDE.md) - How to use COSAP-NX
- [API Reference](docs/API_REFERENCE.md) - Detailed class and method documentation

## Roadmap

- [x] v0.1: BAM → DeepVariant → VCF
- [ ] v0.2: FASTQ → BAM (alignment)
- [ ] v0.3: Multiple callers (HaplotypeCaller)
- [ ] v0.4: Annotation (VEP/ANNOVAR)
