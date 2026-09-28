import pytest

from soysnp_format import formats
from soysnp_format.formats import normalize_genotype, format_genotype


@pytest.mark.parametrize("token,expected", [
    ("AA", ("A", "A")),
    ("AB", ("A", "B")),
    ("BB", ("B", "B")),
    ("AG", ("A", "G")),
    ("A/A", ("A", "A")),
    ("A|B", ("A", "B")),
    ("A B", ("A", "B")),
    ("ab", ("A", "B")),
    ("ct", ("C", "T")),
    ("--", None),
    ("-", None),
    ("", None),
    ("NA", None),
    ("No Call", None),
    ("NN", None),
    (None, None),
])
def test_normalize_genotype(token, expected):
    assert normalize_genotype(token) == expected


@pytest.mark.parametrize("token", ["A", "ABC", "A1", "ZZ", "A/B/C"])
def test_normalize_genotype_rejects(token):
    with pytest.raises(ValueError):
        normalize_genotype(token)


def test_format_genotype_sorts_alleles():
    assert format_genotype("G", "A") == "AG"
    assert format_genotype("B", "A") == "AB"
    assert format_genotype(None, "A") == "--"
    assert format_genotype("A", None, missing="NA") == "NA"


def test_genotype_letters_skips_missing():
    assert formats.genotype_letters([("A", "B"), None, ("B", "B")]) == {"A", "B"}
