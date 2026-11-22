"""
COSAP-NX: COSAP-like pipeline builder with Nextflow backend.

v0.3.0 - COSAP-style step graph, BWA-MEM2, comparison after a run
"""

from .api import BamReader, FastqReader, Trimmer, Mapper, MDUP, Recalibrator, VariantCaller, Pipeline
from .runner import PipelineRunner
from .comparator import VariantComparator

__version__ = "0.3.0"
__all__ = [
    "BamReader",
    "FastqReader",
    "Trimmer",
    "Mapper",
    "MDUP",
    "Recalibrator",
    "VariantCaller",
    "Pipeline",
    "PipelineRunner",
    "VariantComparator",
]
