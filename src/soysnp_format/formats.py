"""Allele-format definitions and genotype token normalization.

GenomeStudio can export the same genotype calls in several allele encodings.
This module covers the two encodings relevant to soybean SNP array work:

AB
    Illumina A/B encoding. Each SNP's two design alleles are labelled A and B
    (the manifest's Allele A / Allele B, defined on the assay design strand),
    and genotypes are reported as AA, AB or BB. The letters A and B are labels,
    not nucleotides.

FORWARD
    Forward-strand nucleotide encoding. Genotypes are reported as the actual
    nucleotides on the forward (+) strand of the reference genome, e.g. AA,
    AG, GG for an A/G SNP. This is the encoding used by most public
    repositories (e.g. SoyBase SNP downloads).

The mapping between the two encodings is SNP-specific: manifest allele A may
correspond to any of the four nucleotides on the forward strand, depending on
the SNP's design strand. Translating therefore requires a per-SNP lookup table
(see :mod:`soysnp_format.lookup`).
"""

AB = "AB"
FORWARD = "FORWARD"

#: Canonical format identifiers accepted on the command line and API.
FORMATS = (AB, FORWARD)

#: Tokens treated as missing / no-call, case-insensitive.
MISSING_TOKENS = frozenset({
    "", "-", "--", "---", "----", "?", ".", "./.",
    "NA", "N/A", "NAN", "NONE", "NULL",
    "NC", "NOCALL", "NO CALL", "NO_CALL", "NOCALL",
    "NN", "N/N", "-/-", "--/--",
})

#: Separators that may appear between the two allele labels of one genotype.
SEPARATORS = ("/", "|", " ", "\t", ":", ";", ",")

#: Nucleotide alphabet used by the FORWARD encoding.
NUCLEOTIDES = frozenset("ACGT")

#: Label alphabet used by the AB encoding.
AB_ALLELES = frozenset("AB")


def normalize_genotype(token):
    """Parse one genotype token into a 2-tuple of allele labels, or None.

    Accepts adjacent labels (``AA``, ``AB``), separated labels (``A/A``,
    ``A B``, ``A|B``), either case, and the missing-call tokens listed in
    :data:`MISSING_TOKENS`. Allele order is preserved as written (callers that
    need a canonical order should sort explicitly).
    """
    if token is None:
        return None
    text = str(token).strip().upper()
    if text in MISSING_TOKENS:
        return None
    for sep in SEPARATORS:
        text = text.replace(sep, "")
    if len(text) != 2:
        raise ValueError(f"cannot parse genotype token {token!r}: expected two allele labels")
    a1, a2 = text[0], text[1]
    if a1 not in "ACGTB" or a2 not in "ACGTB":
        raise ValueError(f"cannot parse genotype token {token!r}: unknown allele labels")
    return (a1, a2)


def genotype_letters(genotypes):
    """Return the set of distinct allele letters observed in an iterable of genotypes.

    Each genotype is a 2-tuple as returned by :func:`normalize_genotype`;
    ``None`` (missing) entries are skipped.
    """
    letters = set()
    for gt in genotypes:
        if gt is None:
            continue
        letters.update(gt)
    return letters


def format_genotype(a1, a2, missing="--"):
    """Render a genotype pair back to a token.

    Alleles are emitted in alphabetical order (``AG``, not ``GA``) so output
    is deterministic regardless of input allele order. ``None`` alleles render
    as the missing token.
    """
    if a1 is None or a2 is None:
        return missing
    return "".join(sorted((a1, a2)))
