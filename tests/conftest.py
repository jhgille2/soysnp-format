"""Shared fixtures for the soysnp-format test suite."""

import pytest

from soysnp_format.lookup import ChipLookup


@pytest.fixture()
def lookup():
    # snp1: A/G SNP, manifest A == forward A  (design strand == forward strand)
    # snp2: T/C SNP, manifest A == forward T
    # snp3: A/T SNP, manifest A == forward A
    # snp4: G/A SNP, manifest A == forward C  (design strand == reverse strand:
    #       manifest A=G -> forward C, manifest B=A -> forward T)
    return ChipLookup(
        {
            "snp1": ("A", "G"),
            "snp2": ("T", "C"),
            "snp3": ("A", "T"),
            "snp4": ("C", "T"),
        },
        name="test-chip",
    )


AB_WIDE = """SNP\tsample1\tsample2\tsample3
snp1\tAA\tAB\tBB
snp2\tAB\tBB\t--
snp3\tAA\tAA\tAB
snp4\tBB\tAB\tAA
"""

FORWARD_WIDE = """SNP\tsample1\tsample2\tsample3
snp1\tAA\tAG\tGG
snp2\tCT\tCC\t--
snp3\tAA\tAA\tAT
snp4\tTT\tCT\tCC
"""

FINAL_REPORT_AB = """[Header]
GSGT Version\t2.0
[Data]
SNP Name\tSample ID\tAllele1 - AB\tAllele2 - AB\tGC Score
snp1\tsample1\tA\tA\t0.9
snp1\tsample2\tA\tB\t0.8
snp2\tsample1\tA\tB\t0.9
snp2\tsample2\tB\tB\t0.7
"""

FINAL_REPORT_FORWARD = """[Header]
GSGT Version\t2.0
[Data]
SNP Name\tSample ID\tAllele1 - Forward\tAllele2 - Forward\tGC Score
snp1\tsample1\tA\tA\t0.9
snp1\tsample2\tA\tG\t0.8
snp2\tsample1\tC\tT\t0.9
snp2\tsample2\tC\tC\t0.7
"""
