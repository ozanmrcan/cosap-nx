"""
COSAP-NX: COSAP-like pipeline builder with Nextflow backend.

v0.1.2 - FASTQ alignment + Germline variant calling (FASTQ/BAM -> DeepVariant -> VCF)
"""

from .api import BamReader, FastqReader, Mapper, VariantCaller, Pipeline
from .runner import PipelineRunner

__version__ = "0.1.2"
__all__ = ["BamReader", "FastqReader", "Mapper", "VariantCaller", "Pipeline", "PipelineRunner"]
