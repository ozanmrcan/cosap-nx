# Changelog

All notable changes to COSAP-NX will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.2] - 2025-12-28

### Added
- FASTQ alignment support using BWA-MEM
- FastqReader class for paired-end FASTQ input configuration
- BWA_ALIGN Nextflow process with automatic sorting and indexing
- Read group handling (@RG tags: ID, SM, PL, LB)
- Mulled container with BWA 0.7.17 + samtools 1.15.1
- germline_from_fastq.py example script
- COSAP capabilities analysis document
- Graduation project roadmap document

### Changed
- VariantCaller now accepts Union[BamReader, List[FastqReader]]
- Updated API documentation with FastqReader usage
- Updated user guide with FASTQ workflow examples
- Resource configuration: BWA process uses 8 CPUs, 12 GB memory

### Technical Details
- Complete FASTQ → BAM → VCF workflow
- BWA-MEM alignment with piped samtools sort
- Backward compatible with v0.1.0 (BAM input still works)
- Branching workflow: FASTQ mode vs BAM mode

### Validated
- Tested with NA12878 chr20 FASTQ data
- Successful alignment and variant calling
- Proper BAM indexing and read group assignment

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
