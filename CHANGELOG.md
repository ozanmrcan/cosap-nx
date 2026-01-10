# Changelog

All notable changes to COSAP-NX will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-01-10

### Added
- **VariantComparator class** - Complete comparison infrastructure matching COSAP API
- **HaplotypeCaller support** - GATK HaplotypeCaller germline variant calling
- **Visual comparison methods**:
  - `draw_venn_diagram()` - Auto-detecting 2-way or 3-way Venn diagrams
  - `draw_venn2_plot()` - 2-way Venn diagram
  - `draw_venn3_plot()` - 3-way Venn diagram
  - `draw_upset_plot()` - UpSet plot for set intersections
  - `draw_similarity_plot()` - Jaccard similarity heatmap
  - `draw_precision_recall_plot()` - Precision/recall scatter plot
- **Metrics computation methods**:
  - `compute_statistics()` - Variant counts (SNPs, indels, totals) per pipeline
  - `compute_overlap()` - Pairwise overlap analysis with Jaccard similarity
  - `compute_metrics_vs_truth()` - Precision, recall, F1 scores vs ground truth
- **Utility methods**:
  - `get_variants()` - Access raw variant DataFrame for custom analysis
  - `create_intersection_bed()` - Export variant sets with set operations (`&`, `|`, `~`)
- **Automatic chromosome name normalization** - Handles chr20 vs 20 naming differences
- Comparison example script (`compare_deepvariant_haplotypecaller.py`)
- Advanced comparison test script (`test_advanced_comparison.py`)

### Changed
- Version bumped to 0.2.0
- Updated dependencies: matplotlib, seaborn, matplotlib-venn, upsetplot, scikit-learn

### Technical Details
- VariantComparator API matches COSAP's VariantComparator exactly
- Uses bcftools for efficient VCF operations (stats, isec)
- Lazy loading of variant data (only when visualization methods called)
- All plots saved as PNG files (non-interactive backend)
- Pipeline naming convention: `{mapper}_{caller}` (e.g., "bwa_deepvariant")

### Validated
- Tested comparison between DeepVariant and HaplotypeCaller on chr20 test data
- Verified all visualization methods generate correct plots
- Confirmed metrics match expected values (precision, recall, F1)
- Validated BED file export with set operations
- Chromosome normalization tested with chr20 vs 20 naming

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
