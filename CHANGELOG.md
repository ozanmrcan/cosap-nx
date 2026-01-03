# Changelog

All notable changes to COSAP-NX will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.2] - 2026-01-04

### Added
- **Mapper class** for alignment-only workflows (COSAP-compatible API)
- **Alignment-only mode**: FASTQ → BAM without variant calling
- **Modular preprocessing architecture**: VariantCaller accepts Mapper input for component chaining
- FASTQ alignment support using BWA-MEM
- FastqReader class for paired-end FASTQ input configuration
- BWA_ALIGN Nextflow process with automatic sorting and indexing
- Read group handling (@RG tags: ID, SM, PL, LB, PU)
- Support for .64 BWA index files (large genomes like hg38)
- Custom read group configuration via Mapper params
- Mulled container with BWA 0.7.17 + samtools 1.15.1
- germline_from_fastq.py example script
- **alignment_only.py example**: Demonstrates FASTQ → BAM workflow
- **mapper_chaining.py example**: Demonstrates Mapper → VariantCaller composition
- COSAP capabilities analysis document
- Graduation project roadmap document

### Changed
- VariantCaller now accepts Union[BamReader, List[FastqReader], Mapper]
- Pipeline.add() now accepts Union[Mapper, VariantCaller]
- Nextflow workflow supports 3 modes: alignment-only, full pipeline, BAM-only
- Updated API documentation with Mapper and modular workflow examples
- Updated user guide with FASTQ workflow examples
- Resource configuration: BWA process uses 8 CPUs, 12 GB memory

### Technical Details
- Complete FASTQ → BAM → VCF workflow
- Alignment-only workflow: FASTQ → BAM (stop before variant calling)
- Component chaining: Mapper output feeds into VariantCaller
- BWA-MEM alignment with piped samtools sort
- Backward compatible with v0.1.0 (BAM input still works)
- Branching workflow: FASTQ mode vs BAM mode vs alignment-only mode
- Read groups fully configurable (SM required, ID/PL/LB/PU optional)

### Validated
- Tested with NA12878 chr20 FASTQ data
- Tested with homework exome data (4.7M read pairs)
- Alignment quality: 98.6% mapping rate, 96.8% properly paired
- Successful alignment and variant calling
- Proper BAM indexing and read group assignment
- Mapper class integration with VariantCaller verified

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
