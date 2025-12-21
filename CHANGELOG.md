# Changelog

All notable changes to COSAP-NX will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2025-12-22

### Added
- Nextflow DSL2 germline variant calling pipeline
- DeepVariant integration (v1.6.0)
- Python API wrapper mimicking COSAP's builder pattern
- BamReader class for BAM input configuration
- VariantCaller class for variant calling configuration
- Pipeline class for workflow orchestration
- PipelineRunner for Nextflow execution
- Comprehensive documentation (SETUP, USER_GUIDE, API_REFERENCE)
- Example germline calling script
- GIAB benchmarking validation (100% SNP recall/precision)

### Technical Details
- Nextflow 25.10.2+ compatibility
- Docker integration for reproducibility
- Configurable CPU/memory resources
- Support for WGS model
- GVCF output option

### Validated Against
- GIAB NA12878 v3.3.2 truth set
- chr20:10MB-10.1MB test region
- hap.py benchmarking tool

### Known Limitations
- Single-sample germline calling only
- DeepVariant only
- No BAM preprocessing (alignment, marking duplicates)
- No joint calling support
- No BED file region filtering
