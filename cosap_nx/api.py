"""
COSAP-like API for building pipelines.

Provides familiar builder classes (BamReader, VariantCaller, Pipeline)
that generate Nextflow-compatible configuration.
"""

import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union


@dataclass
class BamReader:
    """
    Represents an input BAM file.

    Similar to COSAP's BamReader - wraps a BAM path and extracts sample name.
    """

    filename: str
    name: Optional[str] = None

    def __post_init__(self):
        self.filename = os.path.abspath(os.path.normpath(self.filename))

        if self.name is None:
            # Extract name from filename (e.g., "NA12878.bam" -> "NA12878")
            self.name = Path(self.filename).stem

    def get_output(self) -> str:
        """Return the path to the BAM file."""
        return self.filename

    def validate(self) -> None:
        """Check that the BAM file exists."""
        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f"BAM file not found: {self.filename}")

        # Check for index
        bai_path1 = f"{self.filename}.bai"
        bai_path2 = self.filename.replace(".bam", ".bai")
        if not os.path.isfile(bai_path1) and not os.path.isfile(bai_path2):
            raise FileNotFoundError(
                f"BAM index not found. Expected: {bai_path1} or {bai_path2}"
            )


@dataclass
class FastqReader:
    """
    Represents an input FASTQ file for alignment.

    For paired-end sequencing, create two FastqReader instances (read=1, read=2).

    Args:
        filename: Path to FASTQ file (.fastq or .fastq.gz)
        read: Read number (1 or 2 for paired-end)
        name: Sample name (auto-extracted from filename if not provided)
        platform: Sequencing platform (default: "ILLUMINA")
    """

    filename: str
    read: int
    name: Optional[str] = None
    platform: str = "ILLUMINA"

    def __post_init__(self):
        # Normalize path
        self.filename = os.path.abspath(os.path.normpath(self.filename))

        # Validate file exists
        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f"FASTQ file not found: {self.filename}")

        # Auto-extract sample name from filename if not provided
        if self.name is None:
            basename = os.path.basename(self.filename)
            # Remove common suffixes: _1.fastq.gz, _R1.fastq.gz, .R1.fastq, etc.
            self.name = re.sub(r'[._][Rr]?[12]\.f(ast)?q(\.gz)?$', '', basename)

        # Validate read number
        if self.read not in [1, 2]:
            raise ValueError(f"Read number must be 1 or 2, got: {self.read}")

    def get_output(self) -> str:
        """Return the path to the FASTQ file."""
        return self.filename

    def validate(self) -> None:
        """Check that the FASTQ file exists and is readable."""
        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f"FASTQ file not found: {self.filename}")


@dataclass
class VariantCaller:
    """
    Configures variant calling for a sample.

    Supports both BAM and FASTQ inputs:
    - BAM: Provide a BamReader instance
    - FASTQ: Provide a list of 2 FastqReader instances (paired-end reads)

    For v0.1.2, only 'deepvariant' is supported for germline calling.
    """

    library: str
    normal_sample: Optional[Union[BamReader, List[FastqReader]]] = None
    germline: Optional[Union[BamReader, List[FastqReader]]] = None  # Alias for normal_sample
    tumor_sample: Optional[Union[BamReader, List[FastqReader]]] = None
    params: Dict[str, Any] = field(default_factory=dict)
    gvcf: bool = False
    name: Optional[str] = None

    # Supported libraries in v0.1.2
    SUPPORTED_LIBRARIES = ["deepvariant"]

    def __post_init__(self):
        # Normalize library name
        self.library = self.library.lower()

        # Handle germline alias
        if self.germline is not None and self.normal_sample is None:
            self.normal_sample = self.germline

        # Validate library
        if self.library not in self.SUPPORTED_LIBRARIES:
            raise ValueError(
                f"Library '{self.library}' not supported in v0.1.2. "
                f"Supported: {self.SUPPORTED_LIBRARIES}"
            )

        # v0.1.2: Only germline calling
        if self.tumor_sample is not None:
            raise ValueError(
                "Somatic calling (tumor_sample) not supported in v0.1.2"
            )

        if self.normal_sample is None:
            raise ValueError("normal_sample (germline) is required")

        # Validate input type
        if isinstance(self.normal_sample, list):
            # FASTQ input: must have exactly 2 reads (paired-end)
            if len(self.normal_sample) != 2:
                raise ValueError(
                    "FASTQ input requires paired-end reads (exactly 2 FastqReader instances)"
                )
            if not all(isinstance(f, FastqReader) for f in self.normal_sample):
                raise TypeError("List input must contain FastqReader objects")

            # Check read numbers are 1 and 2
            reads = sorted([f.read for f in self.normal_sample])
            if reads != [1, 2]:
                raise ValueError("FASTQ reads must be numbered 1 and 2")

            # Verify both FASTQs have the same sample name
            names = set(f.name for f in self.normal_sample)
            if len(names) > 1:
                raise ValueError(
                    f"FASTQ files have mismatched sample names: {names}. "
                    "Both read 1 and read 2 must have the same 'name' parameter."
                )

        elif not isinstance(self.normal_sample, BamReader):
            raise TypeError(
                "normal_sample must be a BamReader or a list of 2 FastqReader objects"
            )

        # Generate name if not provided
        if self.name is None:
            sample_name = self._get_sample_name()
            self.name = f"{sample_name}_{self.library}"

    def _get_sample_name(self) -> str:
        """Get the sample name from either BAM or FASTQ input."""
        if isinstance(self.normal_sample, BamReader):
            return self.normal_sample.name
        else:  # List[FastqReader]
            # Both should have the same name (validated in __post_init__)
            return self.normal_sample[0].name

    def get_sample_id(self) -> str:
        """Get the sample identifier."""
        return self.params.get(
            "germline_sample_name",
            self._get_sample_name()
        )

    def get_model_type(self) -> str:
        """Get DeepVariant model type (WGS, WES, PACBIO)."""
        return self.params.get("model_type", "WGS")


@dataclass
class Pipeline:
    """
    Collects pipeline steps and builds configuration.

    Usage:
        pipeline = Pipeline()
        pipeline.add(variant_caller)
        config_path = pipeline.build(workdir="/path/to/workdir")
    """

    ref_fasta: Optional[str] = None
    steps: List[VariantCaller] = field(default_factory=list)
    _workdir: Optional[str] = None

    def add(self, step: VariantCaller) -> "Pipeline":
        """
        Add a step to the pipeline.

        Args:
            step: A VariantCaller instance.

        Returns:
            Self for method chaining.
        """
        if not isinstance(step, VariantCaller):
            raise TypeError(
                f"Expected VariantCaller, got {type(step).__name__}"
            )
        self.steps.append(step)
        return self

    def set_reference(self, ref_fasta: str) -> "Pipeline":
        """
        Set the reference genome.

        Args:
            ref_fasta: Path to reference FASTA file.

        Returns:
            Self for method chaining.
        """
        self.ref_fasta = os.path.abspath(os.path.normpath(ref_fasta))
        return self

    def build(
        self,
        workdir: str,
        cpus: int = 8,
        memory: str = "16 GB",
    ) -> str:
        """
        Build the pipeline configuration.

        Generates a params.json file that Nextflow can consume.

        Args:
            workdir: Working directory for outputs.
            cpus: Number of CPUs for DeepVariant.
            memory: Memory allocation (e.g., "16 GB").

        Returns:
            Path to the generated params.json file.

        Raises:
            ValueError: If pipeline is not properly configured.
        """
        if not self.steps:
            raise ValueError("Pipeline has no steps. Use pipeline.add(step).")

        if self.ref_fasta is None:
            raise ValueError(
                "Reference genome not set. Use pipeline.set_reference(path)."
            )

        # Validate reference exists
        if not os.path.isfile(self.ref_fasta):
            raise FileNotFoundError(f"Reference not found: {self.ref_fasta}")

        # For v0.1, we only support a single VariantCaller step
        if len(self.steps) > 1:
            raise ValueError(
                "v0.1 supports only a single VariantCaller step"
            )

        step = self.steps[0]

        # Validate input files
        input_sample = step.normal_sample
        if isinstance(input_sample, BamReader):
            # BAM input - validate BAM file
            input_sample.validate()
        else:  # List[FastqReader]
            # FASTQ input - validate both FASTQ files
            for fastq in input_sample:
                fastq.validate()

        # Setup directories
        self._workdir = os.path.abspath(workdir)
        outdir = os.path.join(self._workdir, "results")
        os.makedirs(outdir, exist_ok=True)

        # Build config based on input type
        config = {
            "sample_id": step.get_sample_id(),
            "ref_fasta": self.ref_fasta,
            "outdir": outdir,
            "model_type": step.get_model_type(),
            "gvcf": step.gvcf,
            "cpus": cpus,
            "memory": memory,
        }

        # Add input-specific parameters
        if isinstance(input_sample, BamReader):
            # BAM input (v0.1 behavior)
            config["bam"] = input_sample.get_output()
        else:  # List[FastqReader]
            # FASTQ input (v0.1.2 behavior)
            fastq_r1 = next(f for f in input_sample if f.read == 1)
            fastq_r2 = next(f for f in input_sample if f.read == 2)
            config["fastq_r1"] = fastq_r1.get_output()
            config["fastq_r2"] = fastq_r2.get_output()

        # Write params.json
        config_path = os.path.join(self._workdir, "params.json")
        with open(config_path, "w") as f:
            json.dump(config, f, indent=2)

        print(f"Pipeline config written to: {config_path}")
        return config_path

    def get_workdir(self) -> Optional[str]:
        """Return the workdir set during build()."""
        return self._workdir
