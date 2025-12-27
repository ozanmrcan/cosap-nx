# COSAP-NX v0.1.2

COSAP-like pipeline builder with **Nextflow** as the execution backend.

## Version History

### v0.1.2 (Current)
**FASTQ → BWA-MEM → sorted BAM → DeepVariant → VCF**

- ✅ FASTQ alignment support (BWA-MEM)
- ✅ Automated BAM sorting and indexing
- ✅ Backward compatible with v0.1.0 (BAM input still works)
- ✅ Comprehensive documentation

### v0.1.0
**BAM → DeepVariant → VCF**

- Minimal vertical slice demonstrating the architecture
- BAM-only input
- DeepVariant germline calling

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

### Option 1: Python API with BAM Input

```python
from cosap_nx import BamReader, VariantCaller, Pipeline, PipelineRunner

# Define BAM input
bam = BamReader("/data/NA12878.bam", name="NA12878")

# Configure variant caller
caller = VariantCaller(
    library="deepvariant",
    normal_sample=bam,
    params={"model_type": "WGS"},
    gvcf=False
)

# Build and run pipeline
pipeline = Pipeline(ref_fasta="/refs/hg38.fa")
pipeline.add(caller)

config_path = pipeline.build(workdir="./workdir")

runner = PipelineRunner()
runner.run_pipeline(config_path)
```

### Option 2: Python API with FASTQ Input (v0.1.2+)

```python
from cosap_nx import FastqReader, VariantCaller, Pipeline, PipelineRunner

# Define FASTQ inputs (paired-end)
fq_r1 = FastqReader("/data/sample_R1.fastq.gz", read=1, name="MySample")
fq_r2 = FastqReader("/data/sample_R2.fastq.gz", read=2, name="MySample")

# Configure variant caller with FASTQ input
caller = VariantCaller(
    library="deepvariant",
    normal_sample=[fq_r1, fq_r2],  # List of FastqReaders
    params={"model_type": "WGS"}
)

# Build and run pipeline
pipeline = Pipeline(ref_fasta="/refs/hg38.fa")  # Must be BWA-indexed!
pipeline.add(caller)

config_path = pipeline.build(workdir="./workdir", cpus=8, memory="16 GB")

runner = PipelineRunner()
runner.run_pipeline(config_path)
```

**Note:** For FASTQ input, the reference must be BWA-indexed:
```bash
bwa index /refs/hg38.fa
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

- [x] v0.1.0: BAM → DeepVariant → VCF
- [x] v0.1.2: FASTQ → BAM (BWA-MEM alignment)
- [ ] v0.2.0: Preprocessing (Mark Duplicates + BQSR)
- [ ] v0.3.0: Multiple callers (HaplotypeCaller, Mutect2)
- [ ] v0.4.0: Annotation (VEP/ANNOVAR)
