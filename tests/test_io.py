from soysnp_format.io import read_genotypes

from conftest import AB_WIDE, FORWARD_WIDE, FINAL_REPORT_AB, FINAL_REPORT_FORWARD


def _write(tmp_path, name, content):
    p = tmp_path / name
    p.write_text(content)
    return str(p)


def test_read_wide_snp_rows(tmp_path):
    data = read_genotypes(_write(tmp_path, "m.tsv", AB_WIDE))
    assert data.layout == "wide"
    assert data.snps == ["snp1", "snp2", "snp3", "snp4"]
    assert data.samples == ["sample1", "sample2", "sample3"]
    assert data.calls["snp1"]["sample1"] == ("A", "A")
    assert data.calls["snp1"]["sample3"] == ("B", "B")
    assert data.calls["snp2"]["sample3"] is None


def test_read_wide_sample_rows(tmp_path):
    # transpose of AB_WIDE: samples as rows, SNPs as columns
    content = "Sample\tsnp1\tsnp2\tsnp3\tsnp4\n" \
              "sample1\tAA\tAB\tAA\tBB\n" \
              "sample2\tAB\tBB\tAA\tAB\n" \
              "sample3\tBB\t--\tAB\tAA\n"
    data = read_genotypes(_write(tmp_path, "m.tsv", content))
    assert data.snps == ["snp1", "snp2", "snp3", "snp4"]
    assert data.samples == ["sample1", "sample2", "sample3"]
    assert data.calls["snp4"]["sample1"] == ("B", "B")


def test_read_wide_orientation_from_chip_ids(tmp_path):
    # header holds the SNP ids; first column holds sample ids
    content = "Sample\tsnp9\tsnp1\tsnp2\n" \
              "sampleA\tAA\tAB\tBB\n" \
              "sampleB\tAB\tAA\tAB\n"
    data = read_genotypes(_write(tmp_path, "m.tsv", content),
                          snp_sets=({"snp1", "snp2", "snp9"},))
    assert data.snps == ["snp9", "snp1", "snp2"]
    assert data.samples == ["sampleA", "sampleB"]
    assert data.calls["snp1"]["sampleA"] == ("A", "B")


def test_read_final_report_ab(tmp_path):
    data = read_genotypes(_write(tmp_path, "r.txt", FINAL_REPORT_AB))
    assert data.layout == "final-report"
    assert data.declared_format == "AB"
    assert data.calls["snp1"]["sample2"] == ("A", "B")
    assert data.calls["snp2"]["sample2"] == ("B", "B")


def test_read_final_report_forward(tmp_path):
    data = read_genotypes(_write(tmp_path, "r.txt", FINAL_REPORT_FORWARD))
    assert data.declared_format == "FORWARD"
    assert data.calls["snp2"]["sample1"] == ("C", "T")


def test_read_final_report_rejects_non_report(tmp_path):
    import pytest
    with pytest.raises(ValueError):
        from soysnp_format.io import read_final_report
        read_final_report(_write(tmp_path, "m.tsv", AB_WIDE))
