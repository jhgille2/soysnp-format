import pytest

from soysnp_format.convert import convert_genotypes, ConversionError, write_wide_matrix
from soysnp_format.io import read_genotypes

from conftest import AB_WIDE, FORWARD_WIDE


def _read(tmp_path, content):
    p = tmp_path / "in.tsv"
    p.write_text(content)
    return read_genotypes(str(p))


def test_ab_to_forward(tmp_path, lookup):
    data = _read(tmp_path, AB_WIDE)
    out = convert_genotypes(data, "AB", "FORWARD", lookup=lookup)
    # snp1: manifest A->A, B->G
    assert out.calls["snp1"]["sample1"] == ("A", "A")
    assert out.calls["snp1"]["sample2"] == ("A", "G")
    assert out.calls["snp1"]["sample3"] == ("G", "G")
    # snp2: manifest A->T, B->C (input allele order preserved internally)
    assert out.calls["snp2"]["sample1"] == ("T", "C")
    assert out.calls["snp2"]["sample2"] == ("C", "C")
    assert out.calls["snp2"]["sample3"] is None
    # snp4: reverse-strand SNP, manifest A->C, B->T
    assert out.calls["snp4"]["sample1"] == ("T", "T")
    assert out.calls["snp4"]["sample2"] == ("C", "T")
    assert out.calls["snp4"]["sample3"] == ("C", "C")


def test_forward_to_ab(tmp_path, lookup):
    data = _read(tmp_path, FORWARD_WIDE)
    out = convert_genotypes(data, "FORWARD", "AB", lookup=lookup)
    assert out.calls["snp1"]["sample2"] == ("A", "B")
    assert out.calls["snp2"]["sample1"] == ("B", "A")
    assert out.calls["snp2"]["sample2"] == ("B", "B")
    assert out.calls["snp2"]["sample3"] is None


def test_round_trip(tmp_path, lookup):
    data = _read(tmp_path, AB_WIDE)
    fwd = convert_genotypes(data, "AB", "FORWARD", lookup=lookup)
    back = convert_genotypes(fwd, "FORWARD", "AB", lookup=lookup)
    for snp in data.snps:
        for sample in data.samples:
            assert back.calls[snp][sample] == data.calls[snp][sample]


def test_auto_detect_source(tmp_path, lookup):
    data = _read(tmp_path, AB_WIDE)
    out = convert_genotypes(data, "auto", "FORWARD", lookup=lookup)
    assert out.calls["snp1"]["sample3"] == ("G", "G")


def test_unknown_snp_error_and_skip(tmp_path, lookup):
    content = "SNP\ts1\nsnp1\tAB\nsnp_UNKNOWN\tAB\n"
    data = _read(tmp_path, content)
    with pytest.raises(ConversionError):
        convert_genotypes(data, "AB", "FORWARD", lookup=lookup)
    out = convert_genotypes(data, "AB", "FORWARD", lookup=lookup, on_unknown_snp="skip")
    assert out.snps == ["snp1"]


def test_bad_allele_error_and_missing(tmp_path, lookup):
    # 'AG' is not valid in AB encoding (enough SNP rows for orientation)
    content = "SNP\ts1\ts2\nsnp1\tAB\tAG\nsnp2\tAA\tBB\nsnp3\tAA\tAA\n"
    data = _read(tmp_path, content)
    with pytest.raises(ConversionError):
        convert_genotypes(data, "AB", "FORWARD", lookup=lookup)
    out = convert_genotypes(data, "AB", "FORWARD", lookup=lookup, on_bad_allele="missing")
    assert out.calls["snp1"]["s2"] is None
    assert out.calls["snp1"]["s1"] == ("A", "G")  # AB -> A/G forward alleles


def test_forward_allele_mismatch_errors(tmp_path, lookup):
    # 'CC' cannot be a forward genotype of snp1 (an A/G SNP)
    content = "SNP\ts1\nsnp1\tCC\n"
    data = _read(tmp_path, content)
    with pytest.raises(ConversionError):
        convert_genotypes(data, "FORWARD", "AB", lookup=lookup)


def test_write_wide_matrix(tmp_path, lookup):
    data = _read(tmp_path, AB_WIDE)
    out = convert_genotypes(data, "AB", "FORWARD", lookup=lookup)
    p = str(tmp_path / "out.tsv")
    write_wide_matrix(out, p)
    lines = open(p).read().splitlines()
    assert lines[0] == "SNP\tsample1\tsample2\tsample3"
    assert lines[1] == "snp1\tAA\tAG\tGG"
