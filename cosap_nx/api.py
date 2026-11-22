"""
COSAP-like API for building pipelines.

Steps point at the step(s) they read from (input_step, germline, tumor), the
same way as in COSAP. Pipeline.build() walks this graph and writes a params.json
where every stage lists its steps by name, and every step names its inputs as
"stage/name". The Nextflow workflow in nf/ runs that file.
"""

import json
import os
import re
from dataclasses import dataclass, field
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import Any, ClassVar, Dict, List, Optional, Union

# Fixed stage order, same as the workflow
STAGES = ["trim", "mapping", "mdup", "calibrate", "variant_calling"]


@dataclass
class FastqReader:
    """
    One FASTQ file. For paired-end reads, create two (read=1 and read=2) with the same name.
    """

    filename: str
    read: int
    name: Optional[str] = None
    platform: str = "ILLUMINA"

    def __post_init__(self):
        self.filename = os.path.abspath(os.path.normpath(self.filename))

        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f"FASTQ file not found: {self.filename}")

        # Sample name from filename if not given: sample_R1.fastq.gz -> sample
        if self.name is None:
            basename = os.path.basename(self.filename)
            self.name = re.sub(r'[._][Rr]?[12]\.f(ast)?q(\.gz)?$', '', basename)

        if self.read not in [1, 2]:
            raise ValueError(f"Read number must be 1 or 2, got: {self.read}")


@dataclass
class BamReader:
    """
    An aligned BAM file (needs a .bai next to it).
    """

    filename: str
    name: Optional[str] = None

    def __post_init__(self):
        self.filename = os.path.abspath(os.path.normpath(self.filename))

        if self.name is None:
            self.name = Path(self.filename).stem

    def validate(self) -> None:
        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f"BAM file not found: {self.filename}")

        bai_path1 = f"{self.filename}.bai"
        bai_path2 = re.sub(r"\.bam$", ".bai", self.filename)
        if not os.path.isfile(bai_path1) and not os.path.isfile(bai_path2):
            raise FileNotFoundError(
                f"BAM index not found. Expected: {bai_path1} or {bai_path2}"
            )


def _reads_name(reads: List[FastqReader]) -> str:
    """Name of a FASTQ pair: the shared sample name, or the names joined."""
    names = list(dict.fromkeys(r.name for r in reads))
    return "-".join(names)


def _check_reads(reads: List[FastqReader]) -> None:
    if not isinstance(reads, list) or not all(isinstance(r, FastqReader) for r in reads):
        raise TypeError("Expected a list of FastqReader")
    if sorted(r.read for r in reads) != [1, 2]:
        raise ValueError("Paired-end input needs one FastqReader with read=1 and one with read=2")


@dataclass
class Trimmer:
    """
    Read trimming (fastp). Input is a FASTQ pair.
    """

    input_step: List[FastqReader]
    name: Optional[str] = None

    STAGE: ClassVar[str] = "trim"

    def __post_init__(self):
        _check_reads(self.input_step)
        if self.name is None:
            self.name = _reads_name(self.input_step)

    def _inputs(self) -> Dict[str, Any]:
        return {"input": self.input_step}

    def _config(self) -> Dict[str, Any]:
        return {}


@dataclass
class Mapper:
    """
    Read alignment. Input is a FASTQ pair or a Trimmer.
    params must have read_groups with at least SM, like in COSAP.
    """

    library: str
    input_step: Union[List[FastqReader], Trimmer]
    params: Dict[str, Any] = field(default_factory=dict)
    name: Optional[str] = None

    STAGE: ClassVar[str] = "mapping"

    def __post_init__(self):
        self.library = self.library.lower()

        if isinstance(self.input_step, Trimmer):
            input_name = self.input_step.name
        else:
            _check_reads(self.input_step)
            input_name = _reads_name(self.input_step)

        if self.name is None:
            self.name = f"{input_name}_{self.library}"

        if "SM" not in self.params.get("read_groups", {}):
            raise ValueError("Please specify a sample name (SM) in params['read_groups'].")

    def _inputs(self) -> Dict[str, Any]:
        return {"input": self.input_step}

    def _config(self) -> Dict[str, Any]:
        return {"library": self.library, "params": self.params}


@dataclass
class MDUP:
    """
    Duplicate marking. spark=True uses GATK MarkDuplicatesSpark, False uses Picard.
    """

    input_step: Union[Mapper, BamReader]
    name: Optional[str] = None
    spark: bool = True

    STAGE: ClassVar[str] = "mdup"

    def __post_init__(self):
        if self.name is None:
            self.name = self.input_step.name

    def _inputs(self) -> Dict[str, Any]:
        return {"input": self.input_step}

    def _config(self) -> Dict[str, Any]:
        return {"spark": self.spark}


@dataclass
class Recalibrator:
    """
    Base quality score recalibration (GATK BQSR).
    """

    input_step: Union[MDUP, Mapper, BamReader]
    name: Optional[str] = None

    STAGE: ClassVar[str] = "calibrate"

    def __post_init__(self):
        if self.name is None:
            self.name = self.input_step.name

    def _inputs(self) -> Dict[str, Any]:
        return {"input": self.input_step}

    def _config(self) -> Dict[str, Any]:
        return {}


@dataclass
class VariantCaller:
    """
    Variant calling. germline is the normal sample, tumor the tumor sample
    (germline callers only use germline). Each is a BAM-producing step or a BamReader.
    """

    library: str
    germline: Optional[Union[Mapper, MDUP, Recalibrator, BamReader]] = None
    tumor: Optional[Union[Mapper, MDUP, Recalibrator, BamReader]] = None
    params: Dict[str, Any] = field(default_factory=dict)
    name: Optional[str] = None
    gvcf: bool = False

    STAGE: ClassVar[str] = "variant_calling"

    def __post_init__(self):
        self.library = self.library.lower()

        if self.germline is None and self.tumor is None:
            raise ValueError("VariantCaller needs a germline and/or tumor input.")

        if self.name is None:
            names = [s.name for s in (self.germline, self.tumor) if s is not None]
            self.name = "_".join(names + [self.library])

    def _inputs(self) -> Dict[str, Any]:
        inputs = {"normal": self.germline, "tumor": self.tumor}
        return {role: step for role, step in inputs.items() if step is not None}

    def _config(self) -> Dict[str, Any]:
        config = {"library": self.library, "params": self.params}
        if self.gvcf:
            config["gvcf"] = True
        return config


Step = Union[Trimmer, Mapper, MDUP, Recalibrator, VariantCaller]


@dataclass
class Pipeline:
    """
    Collects pipeline steps and writes params.json.

    Usage:
        mapper = Mapper(library="bwa", input_step=reads, params={"read_groups": {"SM": "normal"}})
        dv = VariantCaller(library="deepvariant", germline=mapper)
        hc = VariantCaller(library="haplotypecaller", germline=mapper)

        config_path = Pipeline(ref_fasta="ref.fasta").add(mapper).add(dv).add(hc).build(workdir="work")

    Steps that feed an added step are included automatically.
    """

    ref_fasta: Optional[str] = None
    steps: List[Step] = field(default_factory=list)

    def add(self, step: Step) -> "Pipeline":
        if not hasattr(step, "STAGE"):
            raise TypeError(f"Not a pipeline step: {type(step).__name__}")
        self.steps.append(step)
        return self

    def set_reference(self, ref_fasta: str) -> "Pipeline":
        self.ref_fasta = ref_fasta
        return self

    def build(self, workdir: str, cpus: int = 4, memory: str = "8 GB") -> str:
        """
        Write <workdir>/params.json and return its path.

        cpus and memory are the defaults for every step.
        """
        if not self.steps:
            raise ValueError("Pipeline has no steps. Use pipeline.add(step).")

        if self.ref_fasta is None:
            raise ValueError(
                "Reference genome not set. Use Pipeline(ref_fasta=path) or pipeline.set_reference(path)."
            )

        ref_fasta = os.path.abspath(os.path.normpath(self.ref_fasta))
        for path in [ref_fasta, f"{ref_fasta}.fai"]:
            if not os.path.isfile(path):
                raise FileNotFoundError(f"Reference file not found: {path}")

        workdir = os.path.abspath(workdir)
        outdir = os.path.join(workdir, "results")

        inputs: Dict[str, Dict] = {}
        stages: Dict[str, Dict] = {stage: {} for stage in STAGES}

        for step in self.steps:
            self._add_to_config(step, inputs, stages)

        config = {
            "version": version("cosap-nx"),
            "created": datetime.now().strftime(r"%Y-%m-%d--%H-%M-%S"),
            "outdir": outdir,
            "defaults": {"cpus": cpus, "memory": memory},
            "references": {"fasta": ref_fasta},
            "inputs": inputs,
        }
        config.update({stage: steps for stage, steps in stages.items() if steps})

        os.makedirs(outdir, exist_ok=True)
        config_path = os.path.join(workdir, "params.json")
        with open(config_path, "w") as f:
            json.dump(config, f, indent=2)

        print(f"Pipeline config written to: {config_path}")
        return config_path

    def _add_to_config(self, step, inputs: Dict, stages: Dict) -> str:
        """Adds a step (and everything it reads from) to the config. Returns its "stage/name"."""
        if isinstance(step, list):
            _check_reads(step)
            reads = {r.read: r.filename for r in step}
            name = _reads_name(step)
            _put(inputs, name, {"type": "fastq", "r1": reads[1], "r2": reads[2]}, "inputs")
            return f"inputs/{name}"

        if isinstance(step, BamReader):
            step.validate()
            _put(inputs, step.name, {"type": "bam", "bam": step.filename}, "inputs")
            return f"inputs/{step.name}"

        entry = {role: self._add_to_config(s, inputs, stages) for role, s in step._inputs().items()}
        entry.update(step._config())
        _put(stages[step.STAGE], step.name, entry, step.STAGE)
        return f"{step.STAGE}/{step.name}"


def _put(section: Dict, name: str, entry: Dict, stage: str) -> None:
    """Adds an entry, two different steps with the same name in a stage is an error."""
    if name in section and section[name] != entry:
        raise ValueError(
            f"Two different steps are named '{name}' in {stage}. Give one of them another name."
        )
    section[name] = entry
