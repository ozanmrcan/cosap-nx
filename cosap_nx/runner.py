"""
Pipeline runner - executes Nextflow workflows.
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional


class NextflowNotFoundError(Exception):
    """Raised when Nextflow is not installed or not in PATH."""
    pass


class PipelineExecutionError(Exception):
    """Raised when the Nextflow pipeline fails."""
    pass


class PipelineRunner:
    """
    Runs COSAP-NX pipelines using Nextflow as the execution backend.
    """

    def __init__(self, nextflow_path: Optional[str] = None):
        """
        Initialize the runner.

        Args:
            nextflow_path: Path to nextflow executable. If None, searches PATH.
        """
        self.nextflow_path = nextflow_path or self._find_nextflow()
        self._nf_pipeline_dir = Path(__file__).parent.parent / "nf"

    def _find_nextflow(self) -> str:
        """Find Nextflow executable in PATH."""
        nf_path = shutil.which("nextflow")
        if nf_path is None:
            raise NextflowNotFoundError(
                "Nextflow not found in PATH. Please install Nextflow: "
                "https://www.nextflow.io/docs/latest/getstarted.html"
            )
        return nf_path

    def run_pipeline(
        self,
        pipeline_config: str,
        profile: str = "docker",
        resume: bool = True,
        work_dir: Optional[str] = None,
    ) -> Path:
        """
        Execute a pipeline from a config file.

        Args:
            pipeline_config: Path to params.json config file.
            profile: Nextflow profile to use (docker, singularity, etc.).
            resume: Whether to resume from cached results.
            work_dir: Nextflow work directory. Defaults to ./work.

        Returns:
            Path to the output directory.

        Raises:
            FileNotFoundError: If config file doesn't exist.
            PipelineExecutionError: If the pipeline fails.
        """
        config_path = Path(pipeline_config)
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")

        # Load config to get outdir
        with open(config_path) as f:
            config = json.load(f)

        outdir = Path(config.get("outdir", "./results"))
        outdir.mkdir(parents=True, exist_ok=True)

        # Build command
        main_nf = self._nf_pipeline_dir / "main.nf"
        cmd = [
            self.nextflow_path,
            "run",
            str(main_nf),
            "-params-file", str(config_path),
            "-profile", profile,
        ]

        if resume:
            cmd.append("-resume")

        if work_dir:
            cmd.extend(["-work-dir", work_dir])

        # Run Nextflow
        print(f"Running: {' '.join(cmd)}")
        print("-" * 60)

        result = subprocess.run(
            cmd,
            cwd=str(config_path.parent),
            stdout=sys.stdout,
            stderr=sys.stderr,
        )

        if result.returncode != 0:
            raise PipelineExecutionError(
                f"Pipeline failed with exit code {result.returncode}"
            )

        # Return path to VCF output
        vcf_dir = outdir / "vcf" / "deepvariant"
        return vcf_dir


def main():
    """CLI entry point for running pipelines directly."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Run COSAP-NX pipeline with Nextflow"
    )
    parser.add_argument(
        "--config", "-c",
        required=True,
        help="Path to params.json config file"
    )
    parser.add_argument(
        "--profile", "-p",
        default="docker",
        help="Nextflow profile (default: docker)"
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="Disable resume from cache"
    )
    parser.add_argument(
        "--work-dir", "-w",
        help="Nextflow work directory"
    )

    args = parser.parse_args()

    runner = PipelineRunner()
    try:
        output_dir = runner.run_pipeline(
            pipeline_config=args.config,
            profile=args.profile,
            resume=not args.no_resume,
            work_dir=args.work_dir,
        )
        print(f"\nOutput: {output_dir}")
    except (NextflowNotFoundError, PipelineExecutionError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
