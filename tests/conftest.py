from pathlib import Path

import pytest

from cosap_nx import VariantComparator

DATA = Path(__file__).parent / "data"

DEEPVARIANT_VCF = DATA / "deepvariant_chr20.vcf.gz"
HAPLOTYPECALLER_VCF = DATA / "haplotypecaller_chr20.vcf.gz"
TRUTH_VCF = DATA / "truth_chr20_10-10.1mb.vcf.gz"


@pytest.fixture
def comparator():
    return VariantComparator(
        variant_files=[
            {"path": str(DEEPVARIANT_VCF), "caller": "deepvariant", "mapper": "bwa"},
            {"path": str(HAPLOTYPECALLER_VCF), "caller": "haplotypecaller", "mapper": "bwa"},
        ]
    )
