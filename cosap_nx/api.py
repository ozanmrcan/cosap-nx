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
class Mapper:
    """
    Configures read alignment from FASTQ to BAM.

    Matches COSAP's Mapper API for compatibility.

    Args:
        library: Aligner to use (currently only 'bwa' supported)
        input_step: List of 2 FastqReader instances (paired-end reads)
        params: Parameters including read_groups (dict with SM, ID, PL, LB, PU fields)
        name: Optional name for this mapping step

    Example:
        >>> fastq1 = FastqReader("sample_R1.fastq.gz", read=1)
        >>> fastq2 = FastqReader("sample_R2.fastq.gz", read=2)
        >>> mapper = Mapper(
        ...     library="bwa",
        ...     input_step=[fastq1, fastq2],
        ...     params={
        ...         "read_groups": {
        ...             "SM": "NA12878",  # Required: sample name
        ...             "ID": "run1",     # Run ID
        ...             "PL": "ILLUMINA", # Platform
        ...             "LB": "lib1"      # Library
        ...         }
        ...     }
        ... )
    """

    library: str
    input_step: Union[List[FastqReader], FastqReader]
    params: Dict[str, Any] = field(default_factory=dict)
    name: Optional[str] = None
    next_step: Optional["VariantCaller"] = None  # Forward reference

    # Supported aligners in v0.1.2
    SUPPORTED_LIBRARIES = ["bwa"]

    def __post_init__(self):
        # Normalize library name
        self.library = self.library.lower()

        # Validate library
        if self.library not in self.SUPPORTED_LIBRARIES:
            raise ValueError(
                f"Library '{self.library}' not supported. "
                f"Supported aligners: {self.SUPPORTED_LIBRARIES}"
            )

        # Normalize input_step to list
        if isinstance(self.input_step, FastqReader):
            # Single-end (wrap in list)
            self.input_step = [self.input_step]

        # Validate input
        if not isinstance(self.input_step, list):
            raise TypeError("input_step must be a FastqReader or list of FastqReaders")

        if not all(isinstance(f, FastqReader) for f in self.input_step):
            raise TypeError("All items in input_step must be FastqReader instances")

        # For paired-end, validate read numbers
        if len(self.input_step) == 2:
            reads = sorted([f.read for f in self.input_step])
            if reads != [1, 2]:
                raise ValueError("Paired-end reads must be numbered 1 and 2")

            # Verify both FASTQs have the same sample name
            names = set(f.name for f in self.input_step)
            if len(names) > 1:
                raise ValueError(
                    f"FASTQ files have mismatched sample names: {names}. "
                    "Both read 1 and read 2 must have the same 'name' parameter."
                )
        elif len(self.input_step) != 1:
            raise ValueError("input_step must have 1 (single-end) or 2 (paired-end) FASTQ files")

        # Validate read_groups parameter (COSAP requires SM at minimum)
        if "read_groups" not in self.params:
            raise ValueError(
                "params must include 'read_groups' dict. "
                "Example: {'read_groups': {'SM': 'sample_name', 'ID': 'run_id'}}"
            )

        read_groups = self.params["read_groups"]
        if not isinstance(read_groups, dict):
            raise TypeError("params['read_groups'] must be a dict")

        if "SM" not in read_groups:
            raise ValueError(
                "params['read_groups'] must include 'SM' (sample name). "
                "Example: {'read_groups': {'SM': 'NA12878', 'ID': 'run1'}}"
            )

        # Generate name if not provided
        if self.name is None:
            sample_name = self.input_step[0].name
            self.name = f"{sample_name}_{self.library}"

    def get_output(self) -> str:
        """
        Get the output BAM path.

        This is where the aligned BAM will be written.
        """
        # Match COSAP-NX output structure: workdir/results/bam/{library}/{name}.bam
        # The actual path will be constructed by Pipeline.build()
        # For now, return a placeholder that indicates this is an alignment output
        return f"{{workdir}}/results/bam/{self.library}/{self.name}.bam"

    def get_sample_name(self) -> str:
        """Get the sample name from read groups."""
        return self.params["read_groups"]["SM"]

    def get_fastq_inputs(self) -> tuple:
        """Get FASTQ R1 and R2 paths (or just R1 for single-end)."""
        if len(self.input_step) == 2:
            fastq_r1 = next(f for f in self.input_step if f.read == 1)
            fastq_r2 = next(f for f in self.input_step if f.read == 2)
            return (fastq_r1.get_output(), fastq_r2.get_output())
        else:
            return (self.input_step[0].get_output(), None)

    def validate(self) -> None:
        """Validate that input FASTQ files exist."""
        for fastq in self.input_step:
            fastq.validate()


@dataclass
class VariantCaller:
    """
    Configures variant calling for a sample.

    Supports BAM, FASTQ, and Mapper inputs:
    - BAM: Provide a BamReader instance
    - FASTQ: Provide a list of 2 FastqReader instances (paired-end reads)
    - Mapper: Provide a Mapper instance (will use its output BAM)

    For v0.1.2, only 'deepvariant' is supported for germline calling.
    """

    library: str
    normal_sample: Optional[Union[BamReader, List[FastqReader], "Mapper"]] = None
    germline: Optional[Union[BamReader, List[FastqReader], "Mapper"]] = None  # Alias for normal_sample
    tumor_sample: Optional[Union[BamReader, List[FastqReader], "Mapper"]] = None
    params: Dict[str, Any] = field(default_factory=dict)
    gvcf: bool = False
    name: Optional[str] = None

    # Supported libraries in v0.1.3
    SUPPORTED_LIBRARIES = ["deepvariant", "haplotypecaller"]

    def __post_init__(self):
        # Normalize library name
        self.library = self.library.lower()

        # Handle germline alias
        if self.germline is not None and self.normal_sample is None:
            self.normal_sample = self.germline

        # Validate library
        if self.library not in self.SUPPORTED_LIBRARIES:
            raise ValueError(
                f"Library '{self.library}' not supported. "
                f"Supported: {self.SUPPORTED_LIBRARIES}"
            )

        # Library-specific parameter validation
        if self.library == "haplotypecaller":
            # HaplotypeCaller doesn't use model_type parameter
            if "model_type" in self.params:
                raise ValueError("HaplotypeCaller does not use 'model_type' parameter")

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

        elif isinstance(self.normal_sample, Mapper):
            # Mapper input: will use the output BAM from alignment
            # Set the next_step link for pipeline chaining
            self.normal_sample.next_step = self

        elif not isinstance(self.normal_sample, BamReader):
            raise TypeError(
                "normal_sample must be a BamReader, a Mapper, or a list of 2 FastqReader objects"
            )

        # Generate name if not provided
        if self.name is None:
            sample_name = self._get_sample_name()
            self.name = f"{sample_name}_{self.library}"

    def _get_sample_name(self) -> str:
        """Get the sample name from BAM, FASTQ, or Mapper input."""
        if isinstance(self.normal_sample, BamReader):
            return self.normal_sample.name
        elif isinstance(self.normal_sample, Mapper):
            return self.normal_sample.get_sample_name()
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

    Supports Mapper (alignment-only) and VariantCaller (calling) steps.

    Usage:
        # Alignment only:
        pipeline = Pipeline(ref_fasta="ref.fasta")
        pipeline.add(mapper)

        # Alignment + calling:
        pipeline.add(mapper).add(variant_caller)

        # Or just calling (BAM input):
        pipeline.add(variant_caller)

        config_path = pipeline.build(workdir="/path/to/workdir")
    """

    ref_fasta: Optional[str] = None
    steps: List[Union[Mapper, VariantCaller]] = field(default_factory=list)
    _workdir: Optional[str] = None

    def add(self, step: Union[Mapper, VariantCaller]) -> "Pipeline":
        """
        Add a step to the pipeline.

        Args:
            step: A Mapper or VariantCaller instance.

        Returns:
            Self for method chaining.
        """
        if not isinstance(step, (Mapper, VariantCaller)):
            raise TypeError(
                f"Expected Mapper or VariantCaller, got {type(step).__name__}"
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

        Supports:
        - Alignment only: Mapper step
        - Calling only: VariantCaller with BAM input
        - Full pipeline: VariantCaller with FASTQ or Mapper input

        Args:
            workdir: Working directory for outputs.
            cpus: Number of CPUs for processing.
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
                "Reference genome not set. Use Pipeline(ref_fasta=path) or pipeline.set_reference(path)."
            )

        # Validate reference exists
        if not os.path.isfile(self.ref_fasta):
            raise FileNotFoundError(f"Reference not found: {self.ref_fasta}")

        # Setup directories
        self._workdir = os.path.abspath(workdir)
        outdir = os.path.join(self._workdir, "results")
        os.makedirs(outdir, exist_ok=True)

        # Determine pipeline mode based on steps
        if len(self.steps) == 1 and isinstance(self.steps[0], Mapper):
            # Alignment-only mode
            return self._build_alignment_only(outdir, cpus, memory)
        elif len(self.steps) == 1 and isinstance(self.steps[0], VariantCaller):
            # Variant calling mode (with optional integrated alignment)
            return self._build_variant_calling(outdir, cpus, memory)
        else:
            # Multiple steps not yet supported in v0.1.2
            raise ValueError(
                "v0.1.2 supports only: (1) single Mapper step, or (2) single VariantCaller step. "
                "Use VariantCaller(normal_sample=mapper) for chained alignment+calling."
            )

    def _build_alignment_only(self, outdir: str, cpus: int, memory: str) -> str:
        """Build configuration for alignment-only pipeline."""
        mapper = self.steps[0]

        # Validate FASTQ inputs
        mapper.validate()

        # Get FASTQ inputs
        fastq_r1, fastq_r2 = mapper.get_fastq_inputs()

        # Get read group information
        read_groups = mapper.params["read_groups"]
        rg_string = f"@RG\\tID:{read_groups.get('ID', mapper.get_sample_name())}"
        rg_string += f"\\tSM:{read_groups['SM']}"
        if "PL" in read_groups:
            rg_string += f"\\tPL:{read_groups['PL']}"
        if "LB" in read_groups:
            rg_string += f"\\tLB:{read_groups['LB']}"
        if "PU" in read_groups:
            rg_string += f"\\tPU:{read_groups['PU']}"

        # Build configuration
        config = {
            "mode": "alignment_only",
            "sample_id": mapper.get_sample_name(),
            "ref_fasta": self.ref_fasta,
            "outdir": outdir,
            "fastq_r1": fastq_r1,
            "fastq_r2": fastq_r2,
            "read_group": rg_string,
            "cpus": cpus,
            "memory": memory,
        }

        # Write params.json
        config_path = os.path.join(self._workdir, "params.json")
        with open(config_path, "w") as f:
            json.dump(config, f, indent=2)

        print(f"Pipeline config written to: {config_path}")
        print(f"Mode: Alignment only (FASTQ → BAM)")
        return config_path

    def _build_variant_calling(self, outdir: str, cpus: int, memory: str) -> str:
        """Build configuration for variant calling pipeline."""
        step = self.steps[0]

        # Validate input files based on input type
        input_sample = step.normal_sample
        if isinstance(input_sample, BamReader):
            input_sample.validate()
        elif isinstance(input_sample, Mapper):
            input_sample.validate()
        else:  # List[FastqReader]
            for fastq in input_sample:
                fastq.validate()

        # Build config based on input type
        config = {
            "sample_id": step.get_sample_id(),
            "ref_fasta": self.ref_fasta,
            "outdir": outdir,
            "variant_caller_library": step.library,
            "gvcf": step.gvcf,
            "cpus": cpus,
            "memory": memory,
        }

        # Add model_type only for DeepVariant
        if step.library == "deepvariant":
            config["model_type"] = step.get_model_type()

        # Add input-specific parameters
        if isinstance(input_sample, BamReader):
            # BAM input (v0.1 behavior)
            config["bam"] = input_sample.get_output()
        elif isinstance(input_sample, Mapper):
            # Mapper input (chained alignment + calling)
            fastq_r1, fastq_r2 = input_sample.get_fastq_inputs()
            config["fastq_r1"] = fastq_r1
            config["fastq_r2"] = fastq_r2
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
