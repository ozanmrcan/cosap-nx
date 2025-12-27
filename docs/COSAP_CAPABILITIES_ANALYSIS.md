# COSAP Comprehensive Capabilities Analysis

**Purpose:** Complete feature inventory of COSAP to guide COSAP-NX development roadmap.

**Source:** COSAP official repository analysis (December 2025)

---

## Table of Contents
1. [Preprocessing Tools](#1-preprocessing-tools--steps)
2. [Variant Calling Tools](#2-variant-calling-tools)
3. [Post-Processing & Analysis](#3-post-processing--analysis)
4. [API Structure](#4-api-structure)
5. [Workflow Features](#5-workflow-features)
6. [GPU/Acceleration Support](#6-gpuacceleration-support)
7. [File Format Support](#7-file-format-support)
8. [COSAP-NX Roadmap](#cosap-nx-development-roadmap)

---

## 1. PREPROCESSING TOOLS & STEPS

### Trimming
| Tool | Implementation | Features |
|------|----------------|----------|
| **Fastp** | `_trimmer.py` | Paired-end trimming, multi-threaded, JSON reports, quality filtering |

### Alignment/Mapping
| Tool | Implementation | Features | Acceleration |
|------|----------------|----------|--------------|
| **BWA** | `_bwa_mapper.py` | BWA-MEM, read groups (ID/SM/LB/PL/PU) | CPU, GPU (Parabricks fq2bam) |
| **BWA-MEM2** | `_bwa2_mapper.py` | Faster BWA-MEM, same read groups | CPU only |
| **Bowtie2** | `_bowtie_mapper.py` | Alternative aligner, read groups | CPU only |

### Sorting
| Tool | Implementation | Features |
|------|----------------|----------|
| **Samtools** | `_sorter.py` | Coordinate/queryname sorting, multi-threaded, piping |

### Duplicate Marking
| Tool | Implementation | Features | Acceleration |
|------|----------------|----------|--------------|
| **Picard MarkDuplicates** | `_mark_duplicate.py` | Mark or delete, metrics file | CPU |
| **GATK MarkDuplicatesSpark** | `_mark_duplicate.py` | Spark parallel processing | Spark |
| **Parabricks MarkDup** | Parabricks | Queryname sort input required | GPU |
| **Elprep** | `_elprep_preprocess.py` | Combined with BQSR, single-pass | CPU (optimized) |

### Base Quality Score Recalibration (BQSR)
| Tool | Implementation | Features | Acceleration |
|------|----------------|----------|--------------|
| **GATK BaseRecalibrator + ApplyBQSR** | `_base_recalibrator.py` | Known sites (Mills, dbSNP, 1000G), BED support | CPU |
| **Parabricks BQSR** | Parabricks | Combined bqsr + applybqsr | GPU |
| **Elprep BQSR** | `_elprep_preprocess.py` | Integrated with markdup | CPU (optimized) |

### Other Preprocessing
- BAM Indexing (Picard BuildBamIndex)
- BAM Merging (Picard MergeSamFiles)
- Split BAM by Chromosome (scatter-gather)

---

## 2. VARIANT CALLING TOOLS

### Germline Callers

| Caller | Implementation | Features |
|--------|----------------|----------|
| **HaplotypeCaller** | `_haplotypecaller_variantcaller.py` | VCF/gVCF, CNN scoring, scatter-gather, BED support, SNP/INDEL separation |
| **DeepVariant** | `_deepvariant_variantcaller.py` | WGS model, VCF/gVCF, multi-threaded sharding, Docker-based |
| **VarScan Germline** | `_varscan_germline_variantcaller.py` | Mpileup-based, configurable p-value |
| **Octopus** | `_octopus_variantcaller.py` | Bayesian haplotype-based, tumor-normal capable |

### Somatic Callers

| Caller | Implementation | Features | GPU Support |
|--------|----------------|----------|-------------|
| **Mutect2** | `_mutect2_variantcaller.py` | Tumor-only or paired, scatter-gather, FilterMutectCalls | ✅ (Parabricks) |
| **Strelka2** | `_strelka_variantcaller.py` | Tumor-normal, separate SNV/INDEL outputs, workflow-based | ❌ |
| **Manta** | `_manta_variantcaller.py` | Structural variants, tumor-normal, BED support | ❌ |
| **VarDict** | `_vardict_variantcaller.py` | Tumor-normal, BED required, AF thresholds | ❌ |
| **VarScan** | `_varscan_variantcaller.py` | Mpileup-based, separate SNV/INDEL, processSomatic filtering | ❌ |
| **MuSE** | `_muse_variantcaller.py` | Two-step (call + sump), dbSNP integration | ❌ |
| **SomaticSniper** | `_somaticsniper_variantcaller.py` | Tumor-normal, LOH detection | ❌ |
| **VarNet** | `_varnet_variantcaller.py` | Deep learning-based, separate SNV/INDEL, Docker-based | ❌ |

**Total:** 4 germline + 8 somatic + 1 SV caller = **13 variant callers**

---

## 3. POST-PROCESSING & ANALYSIS

### VCF Filtering
- GATK FilterVariantTranches (HaplotypeCaller)
- GATK FilterMutectCalls (Mutect2)
- GATK SelectVariants (SNP/INDEL extraction)
- VarScan processSomatic

### Annotation Tools

| Annotator | Implementation | Features | Databases |
|-----------|----------------|----------|-----------|
| **Ensembl VEP** | `_ensembl_vep_annotator.py` | Offline cache, multi-threaded, Docker | Ensembl |
| **ANNOVAR** | `_annovar_annotator.py` | Multi-protocol annotation | ensGene, gnomAD 4.0, ICGC 28, avSNP150, dbNSFP 4.2a, ClinVar |
| **InterVar** | `_intervar_annotator.py` | ACMG interpretation, pathogenicity | InterVar DB |
| **CancerVar** | `_cancervar_annotatator.py` | Cancer-specific, somatic classification | COSMIC |
| **PharmCAT** | `_pharmcat_annotator.py` | Pharmacogenomics, JSON output | PharmGKB |
| **AnnotSV** | `_annotsv_annotator.py` | Structural variant annotation, TSV output | AnnotSV DB |

### Quality Control

| Tool | Implementation | Features |
|------|----------------|----------|
| **Qualimap** | `_qualimap.py` | BAM QC, HTML reports, BED support, multi-threaded |
| **Mosdepth** | `_mosdepth.py` | Fast coverage analysis, per-base/region depth, multi-threaded |

### Specialized Analyses

| Analysis Type | Tool | Implementation | Features |
|---------------|------|----------------|----------|
| **CNV Calling** | CNVkit | `_cnvkit_cnv_caller.py` | Tumor-normal, gene-level metrics, BED targets |
| **MSI Detection** | MSIsensor-pro | `_msisensorpro_msicaller.py` | Tumor-normal, pre-computed microsatellite sites |
| **Gene Fusion** | GeneFuse | `_genefuse_fusion_caller.py` | Direct from FASTQ, JSON/HTML output, multi-threaded |

---

## 4. API STRUCTURE

### Builder Classes (Low-Level API)

**File Readers:**
- `FastqReader` - Paired-end FASTQ files
- `BamReader` - BAM files
- `VCFReader` - VCF files

**Preprocessing:**
- `Trimmer` - Read trimming configuration
- `Mapper` - Alignment (BWA/BWA2/Bowtie2)
- `Sorter` - BAM sorting
- `MDUP` - Duplicate marking
- `Recalibrator` - BQSR configuration
- `Elprep` - Combined preprocessing
- `Indexer` - BAM indexing
- `Merger` - BAM merging

**Variant Calling:**
- `VariantCaller` - Unified caller interface
  - Germline/tumor inputs
  - VCF/gVCF modes
  - BED file intervals
  - Library selection

**Analysis:**
- `Annotator` - Variant annotation
- `QualityController` - QC analysis
- `CNVCaller` - CNV detection
- `MSICaller` - MSI detection
- `GeneFusionCaller` - Fusion detection

### High-Level APIs

**DNAPipeline:**
- `DNAPipelineInput` dataclass:
  ```python
  ANALYSIS_TYPE: "somatic" | "germline"
  WORKDIR: str
  NORMAL_SAMPLE: tuple[str, str]  # (R1, R2)
  TUMOR_SAMPLES: List[Tuple[str, str]]
  BED_FILE: str
  MAPPERS: List[str]
  VARIANT_CALLERS: List[str]
  ANNOTATORS: List[str]
  NORMAL_SAMPLE_NAME: str
  TUMOR_SAMPLE_NAME: str
  BAM_QC: str  # "qualimap" or "mosdepth"
  GVCF: bool
  MSI: bool
  GENEFUSION: bool
  DEVICE: str  # "cpu" or "gpu"
  ```
- Automatic pipeline construction
- Multi-mapper, multi-caller combinations
- Built-in preprocessing workflow

**VariantMultipleAnnotator:**
- Multi-tool annotation (VEP + InterVar + CancerVar)
- Pandas DataFrame output
- Result merging

---

## 5. WORKFLOW FEATURES

### Multi-Sample Support
- ✅ Multiple tumor samples per pipeline
- ✅ Tumor-normal pairing
- ✅ Sample name tracking via read groups
- ✅ Separate output per sample

### Joint Calling
- ❌ **Not implemented** (no GenotypeGVCFs or CombineGVCFs)
- ⚠️ gVCF generation supported (HaplotypeCaller, DeepVariant)
- ⚠️ gVCF can be used for downstream joint calling externally

### Parallelization Strategies

**Scatter-Gather** (`_scatter_gather.py`):
- Interval-based splitting (BED file or auto-generated)
- Parallel variant calling per interval
- VCF/gVCF merging (GATK MergeVcfs, SortVcf)
- Thread-aware configuration

**Multi-Threading:**
- Per-tool thread configuration
- ProcessPoolExecutor for parallel execution
- Spark-based parallelization (MarkDuplicatesSpark)

**Snakemake Orchestration:**
- Rules in `cosap/snakemake_workflows/rules/`
- Dependency resolution
- Checkpoint-based execution

---

## 6. GPU/ACCELERATION SUPPORT

### Parabricks (NVIDIA Clara)
**Docker Image:** `nvcr.io/nvidia/clara/clara-parabricks:4.2.0-1`

| Tool | Command | Function |
|------|---------|----------|
| **fq2bam** | `pbrun fq2bam` | BWA-MEM alignment (GPU) |
| **markdup** | `pbrun markdup` | Duplicate marking (GPU) |
| **bamsort** | `pbrun bamsort` | BAM sorting (GPU) |
| **bqsr** | `pbrun bqsr` | BQSR (GPU) |
| **applybqsr** | `pbrun applybqsr` | Apply BQSR (GPU) |
| **mutectcaller** | `pbrun mutectcaller` | Mutect2 (GPU) |

### Elprep
- CPU-based high-performance preprocessing
- Single-pass BAM processing
- Combined duplicate marking + BQSR

### GATK Spark
- MarkDuplicatesSpark
- Local Spark master mode
- No distributed cluster support

---

## 7. FILE FORMAT SUPPORT

### Input Formats
| Format | Support Level | Usage |
|--------|---------------|-------|
| **FASTQ** | ✅ Full | Paired-end, gzip-compressed |
| **BAM** | ✅ Full | Aligned reads, indexed |
| **CRAM** | ⚠️ Implicit | Via Samtools |
| **VCF** | ✅ Full | Variant input for annotation |
| **BED** | ✅ Full | Interval/target regions |
| **AVinput** | ✅ Full | ANNOVAR format |

### Output Formats
| Format | Tools | Usage |
|--------|-------|-------|
| **VCF** | All callers | Standard variant calls |
| **gVCF** | HaplotypeCaller, DeepVariant | Genomic VCF |
| **BAM** | Preprocessing | Aligned/processed reads |
| **JSON** | GeneFuse, PharmCAT | Gene fusion, pharmacogenomics |
| **TSV** | VEP, AnnotSV | Annotation tables |
| **TXT** | ANNOVAR, InterVar, CancerVar, MSI, CNV | Text reports |
| **HTML** | Qualimap, GeneFuse | Quality reports |

### Reference Genome
- **Primary:** hg38
- **Formats:** FASTA, BWA index, Bowtie2 index
- **Known Sites:** Mills indels, dbSNP, 1000 Genomes

---

## SUMMARY STATISTICS

### Tool Inventory
- **3 Mappers:** BWA, BWA-MEM2, Bowtie2
- **13 Variant Callers:** 4 germline, 8 somatic, 1 SV
- **6 Annotators:** VEP, ANNOVAR, InterVar, CancerVar, PharmCAT, AnnotSV
- **2 QC Tools:** Qualimap, Mosdepth
- **1 Trimmer:** Fastp
- **4 Duplicate Markers:** Picard, GATK Spark, Parabricks, Elprep
- **1 CNV Caller:** CNVkit
- **1 MSI Caller:** MSIsensor-pro
- **1 Gene Fusion Caller:** GeneFuse

### Key Strengths
✅ Comprehensive somatic variant calling suite
✅ GPU acceleration support (Parabricks)
✅ Flexible API with builder pattern
✅ Multi-tool comparison capability
✅ Scatter-gather parallelization
✅ Docker-based tool isolation
✅ Snakemake workflow orchestration
✅ Combined preprocessing optimization (Elprep)

### Notable Limitations
❌ No joint calling implementation (GenotypeGVCFs)
❌ No SnpEff annotator
❌ GPU support limited to Parabricks tools
❌ No distributed computing (only local Spark)
⚠️ hg38 reference only (no hg19/GRCh37 detected)

---

## COSAP-NX Development Roadmap

### ✅ Completed (v0.1-v0.2)
- [x] BAM → DeepVariant → VCF (v0.1)
- [x] FASTQ → BWA-MEM → sorted BAM → DeepVariant → VCF (v0.2)
- [x] Python API (BamReader, FastqReader, VariantCaller, Pipeline)
- [x] Nextflow execution engine
- [x] Docker containerization
- [x] Benchmarking against GIAB (100% SNP recall/precision)

### 🎯 Short-Term Goals (v0.3-v0.4)

**v0.3: Complete Preprocessing**
- [ ] Trimming (Fastp)
- [ ] Mark Duplicates (Picard)
- [ ] BQSR (GATK BaseRecalibrator + ApplyBQSR)
- [ ] Full BAM preprocessing pipeline: FASTQ → trimmed → aligned → sorted → dedup → BQSR → BAM

**v0.4: Multiple Variant Callers**
- [ ] HaplotypeCaller (germline)
- [ ] Mutect2 (somatic)
- [ ] Strelka2 (somatic)
- [ ] Multi-caller support in VariantCaller builder

### 🔮 Medium-Term Goals (v0.5-v0.7)

**v0.5: Annotation**
- [ ] Ensembl VEP integration
- [ ] ANNOVAR integration
- [ ] Annotator builder class
- [ ] Multi-annotator support

**v0.6: Quality Control & Analysis**
- [ ] Qualimap integration
- [ ] Mosdepth integration
- [ ] QualityController builder
- [ ] CNVkit integration
- [ ] MSIsensor-pro integration

**v0.7: Advanced Features**
- [ ] Scatter-gather parallelization
- [ ] BED file interval support
- [ ] gVCF generation
- [ ] Tumor-normal pairing

### 🚀 Long-Term Goals (v1.0+)

**v1.0: COSAP Feature Parity**
- [ ] All 13 variant callers
- [ ] All 6 annotators
- [ ] All preprocessing options
- [ ] Multi-sample workflows
- [ ] Somatic calling workflows
- [ ] Gene fusion detection

**v2.0: Advanced Acceleration**
- [ ] GPU support (Parabricks integration?)
- [ ] Elprep integration
- [ ] AWS Batch / cloud execution
- [ ] Distributed computing support

### 📋 Priority Matrix

| Feature | Priority | Complexity | Impact |
|---------|----------|------------|--------|
| Mark Duplicates | **HIGH** | Low | High |
| BQSR | **HIGH** | Medium | High |
| HaplotypeCaller | **HIGH** | Medium | High |
| Trimming | Medium | Low | Medium |
| Mutect2 | Medium | Medium | High |
| VEP Annotation | Medium | High | Medium |
| Scatter-gather | Low | High | High |
| GPU Support | Low | Very High | Medium |

---

**Document Version:** 1.0
**Last Updated:** December 2025
**Author:** COSAP-NX Development Team
**Source Repository:** `/home/ozan/grad_project/cosap_official_repo`
