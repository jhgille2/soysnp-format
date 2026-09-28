"""Tests for the bundled SoySNP50K / BARCSoySNP6K lookup tables.

These pin the derivation results: row counts, the Illumina A/B rule spot
checks, ss-ID aliasing, and the empirical chip property that no A/T or C/G
SNPs occur (so every mapping is unambiguous).
"""

import pytest

from soysnp_format.convert import convert_genotypes
from soysnp_format.io import read_genotypes
from soysnp_format.lookup import BUNDLED_CHIPS, load_lookup


@pytest.fixture(scope="module")
def lookup50k():
    return load_lookup("soysnp50k")


@pytest.fixture(scope="module")
def lookup6k():
    return load_lookup("soysnp6k")


def test_bundled_chip_names():
    assert set(BUNDLED_CHIPS) == {"soysnp50k", "soysnp6k"}


def test_row_counts(lookup50k, lookup6k):
    assert len(lookup50k) == 60800
    assert len(lookup6k) == 5989


def test_spot_check_ag_snp(lookup50k):
    # BARC_1.01_Gm01_2033_G_A: A/G SNP -> Illumina allele A = A
    assert lookup50k.alleles("BARC_1.01_Gm01_2033_G_A") == ("A", "G")


def test_spot_check_tc_snp(lookup50k):
    # BARC_1.01_Gm01_25990_C_T: C/T SNP -> Illumina allele A = T (Table 1 rule)
    assert lookup50k.alleles("BARC_1.01_Gm01_25990_C_T") == ("T", "C")


def test_ss_id_alias(lookup50k):
    assert "ss715578672" in lookup50k
    assert lookup50k.alleles("ss715578672") == lookup50k.alleles("BARC_1.01_Gm01_2033_G_A")


def test_no_ambiguous_snp_types(lookup50k):
    # Empirical chip property: SoySNP50K contains no A/T or C/G SNPs, so the
    # Illumina A/B rule needs no sequence walking anywhere on the chip.
    bad = [s for s in lookup50k.mapping
           if set(lookup50k.alleles(s)) in ({"A", "T"}, {"C", "G"})]
    assert bad == []


def test_6k_subset_of_50k(lookup50k, lookup6k):
    for snp, pair in lookup6k.mapping.items():
        assert lookup50k.alleles(snp) == pair


def test_convert_with_bundled_lookup(tmp_path, lookup50k):
    content = ("SNP\ts1\ts2\n"
               "BARC_1.01_Gm01_2033_G_A\tAA\tAB\n"
               "BARC_1.01_Gm01_25990_C_T\tAB\tBB\n")
    p = tmp_path / "in.tsv"
    p.write_text(content)
    data = read_genotypes(str(p))
    out = convert_genotypes(data, "auto", "FORWARD", lookup=lookup50k)
    assert out.calls["BARC_1.01_Gm01_2033_G_A"]["s2"] == ("A", "G")
    # C/T SNP: AB-format A == forward T, AB-format B == forward C
    assert out.calls["BARC_1.01_Gm01_25990_C_T"]["s1"] == ("T", "C")
    assert out.calls["BARC_1.01_Gm01_25990_C_T"]["s2"] == ("C", "C")
    back = convert_genotypes(out, "FORWARD", "AB", lookup=lookup50k)
    assert back.calls["BARC_1.01_Gm01_2033_G_A"]["s2"] == ("A", "B")
