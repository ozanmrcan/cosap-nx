"""
COSAP-NX: COSAP-like pipeline builder with Nextflow backend.

v0.1 - Germline variant calling (BAM -> DeepVariant -> VCF)
"""

from .api import BamReader, VariantCaller, Pipeline
from .runner import PipelineRunner

__version__ = "0.1.0"
__all__ = ["BamReader", "VariantCaller", "Pipeline", "PipelineRunner"]
