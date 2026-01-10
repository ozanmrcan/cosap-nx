# COSAP-NX v0.2.0

COSAP-like pipeline builder with **Nextflow** as the execution backend.

## Version History

### v0.2.0 (Current)
**Multi-caller comparison + HaplotypeCaller support**

- ✅ **VariantComparator** - Complete comparison infrastructure matching COSAP API
- ✅ **HaplotypeCaller** - GATK HaplotypeCaller germline variant calling
- ✅ **Visual comparison** - Venn diagrams, UpSet plots, similarity heatmaps, precision/recall plots
- ✅ **Metrics computation** - Variant statistics, overlap analysis, truth comparison
- ✅ **BED export** - Set operations for variant filtering
- ✅ Automatic chromosome name normalization

### v0.1.2
**FASTQ → BWA-MEM → sorted BAM → DeepVariant → VCF**

- ✅ FASTQ alignment support (BWA-MEM)
- ✅ Automated BAM sorting and indexing
- ✅ Backward compatible with v0.1.0 (BAM input still works)

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

### Option 1: Compare Multiple Variant Callers

```python
from cosap_nx import BamReader, VariantCaller, Pipeline, PipelineRunner, VariantComparator

# Run DeepVariant pipeline
bam = BamReader("/data/NA12878.bam", name="NA12878")
caller_dv = VariantCaller(library="deepvariant", normal_sample=bam)
pipeline_dv = Pipeline(ref_fasta="/refs/hg38.fa")
pipeline_dv.add(caller_dv)
config_dv = pipeline_dv.build(workdir="./workdir_dv")
PipelineRunner().run_pipeline(config_dv)

# Run HaplotypeCaller pipeline
caller_hc = VariantCaller(library="haplotypecaller", normal_sample=bam)
pipeline_hc = Pipeline(ref_fasta="/refs/hg38.fa")
pipeline_hc.add(caller_hc)
config_hc = pipeline_hc.build(workdir="./workdir_hc")
PipelineRunner().run_pipeline(config_hc)

# Compare results
comparator = VariantComparator(
    variant_files=[
        {"path": "./workdir_dv/results/vcf/deepvariant/NA12878.vcf.gz",
         "caller": "deepvariant", "mapper": "bwa"},
        {"path": "./workdir_hc/results/vcf/haplotypecaller/NA12878.vcf.gz",
         "caller": "haplotypecaller", "mapper": "bwa"},
    ]
)

# Generate comparison metrics and visualizations
comparator.compute_statistics(output_dir="./comparison_results")
comparator.compute_overlap(output_dir="./comparison_results")
comparator.compute_metrics_vs_truth(
    truth_vcf="/data/truth.vcf.gz",
    output_dir="./comparison_results"
)
comparator.draw_venn_diagram(output_file="./comparison_results/venn.png")
comparator.draw_upset_plot(output_file="./comparison_results/upset.png")
comparator.draw_similarity_plot(output_file="./comparison_results/similarity.png")
comparator.draw_precision_recall_plot(
    truth_vcf="/data/truth.vcf.gz",
    output_file="./comparison_results/pr.png"
)
```

### Option 2: Python API with BAM Input

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

### Option 3: Python API with FASTQ Input

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

### Option 4: Direct Nextflow (for testing)

```bash
nextflow run nf/main.nf \
    -params-file examples/sample_params.json \
    -profile docker
```

### Option 5: CLI Runner

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

## Variant Callers

### DeepVariant
- **Library**: `"deepvariant"`
- **Model types**: WGS, WES, PACBIO
- **Parameters**: `model_type`, `gvcf`

### HaplotypeCaller
- **Library**: `"haplotypecaller"`
- **Parameters**: `gvcf`, `germline_sample_name`

## Comparison Features

COSAP-NX provides complete variant comparison matching COSAP's API:

### Metrics
- `compute_statistics()` - Variant counts (SNPs, indels, totals)
- `compute_overlap()` - Pairwise overlap and Jaccard similarity
- `compute_metrics_vs_truth()` - Precision, recall, F1 scores

### Visualizations
- `draw_venn_diagram()` - Auto 2-way or 3-way Venn diagrams
- `draw_upset_plot()` - UpSet plot for set intersections
- `draw_similarity_plot()` - Jaccard similarity heatmap
- `draw_precision_recall_plot()` - Performance comparison scatter plot

### Utilities
- `get_variants()` - Access raw variant DataFrame
- `create_intersection_bed()` - Export variant sets with set operations (`&`, `|`, `~`)

See `examples/compare_deepvariant_haplotypecaller.py` for complete comparison workflow.

## Configuration

### Variant Caller Parameters

| Parameter | DeepVariant | HaplotypeCaller | Description |
|-----------|-------------|-----------------|-------------|
| `model_type` | ✅ | ❌ | WGS, WES, or PACBIO |
| `gvcf` | ✅ | ✅ | Output gVCF format |
| `germline_sample_name` | ❌ | ✅ | Sample name in VCF |
| `cpus` | ✅ | ✅ | CPU cores (default: 8) |
| `memory` | ✅ | ✅ | Memory allocation (default: 16 GB) |

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
- [x] v0.2.0: Multi-caller comparison + HaplotypeCaller
- [ ] v0.3.0: Additional variant callers (Mutect2, Strelka)
- [ ] v0.4.0: Preprocessing (Mark Duplicates + BQSR)
- [ ] v0.5.0: Annotation (VEP/ANNOVAR)
