from soysnp_format.detect import detect_format, AB, FORWARD, MIXED, INDETERMINATE
from soysnp_format.io import read_genotypes

from conftest import AB_WIDE, FORWARD_WIDE, FINAL_REPORT_AB, FINAL_REPORT_FORWARD


def _read(tmp_path, content, **kw):
    p = tmp_path / "in.tsv"
    p.write_text(content)
    return read_genotypes(str(p), **kw)


def test_detect_ab_wide(tmp_path, lookup):
    result = detect_format(_read(tmp_path, AB_WIDE), lookup=lookup)
    assert result.verdict == AB
    assert result.evidence["n_ab_evidence"] == 4
    assert result.evidence["n_forward_evidence"] == 0


def test_detect_forward_wide(tmp_path, lookup):
    result = detect_format(_read(tmp_path, FORWARD_WIDE), lookup=lookup)
    assert result.verdict == FORWARD
    assert result.evidence["n_forward_evidence"] == 4


def test_detect_ab_final_report(tmp_path, lookup):
    result = detect_format(_read(tmp_path, FINAL_REPORT_AB), lookup=lookup)
    assert result.verdict == AB
    assert result.declared_format == AB


def test_detect_forward_final_report(tmp_path, lookup):
    result = detect_format(_read(tmp_path, FINAL_REPORT_FORWARD), lookup=lookup)
    assert result.verdict == FORWARD


def test_detect_declared_header_disagreement_noted(tmp_path, lookup):
    # data are forward nucleotides but headers claim AB
    content = FINAL_REPORT_FORWARD.replace("Allele1 - Forward", "Allele1 - AB") \
                                  .replace("Allele2 - Forward", "Allele2 - AB")
    result = detect_format(_read(tmp_path, content), lookup=lookup)
    assert result.verdict == FORWARD
    assert any("headers declare" in n.lower() for n in result.notes)


def test_detect_mixed(tmp_path, lookup):
    content = ("SNP\ts1\ts2\n"
               "snp1\tAB\tAA\n"     # B call -> AB evidence
               "snp2\tCT\tCC\n")    # C/T calls -> forward evidence
    result = detect_format(_read(tmp_path, content), lookup=lookup)
    assert result.verdict == MIXED


def test_detect_indeterminate_all_A(tmp_path, lookup):
    content = "SNP\ts1\ts2\nsnp1\tAA\tAA\nsnp3\tAA\t--\n"
    result = detect_format(_read(tmp_path, content), lookup=lookup)
    assert result.verdict == INDETERMINATE


def test_detect_without_lookup(tmp_path):
    result = detect_format(_read(tmp_path, AB_WIDE), lookup=None, chip="none")
    assert result.verdict == AB
