"""
COSAP-NX: COSAP-like pipeline builder with Nextflow backend.

v0.2.0 - Comparison infrastructure + HaplotypeCaller support
"""

from .api import BamReader, FastqReader, Mapper, VariantCaller, Pipeline
from .runner import PipelineRunner
from .comparator import VariantComparator

__version__ = "0.2.0"
__all__ = [
    "BamReader",
    "FastqReader",
    "Mapper",
    "VariantCaller",
    "Pipeline",
    "PipelineRunner",
    "VariantComparator",
]
