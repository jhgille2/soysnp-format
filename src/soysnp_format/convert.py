"""Translate genotypes between A/B and FORWARD encodings.

Translation is SNP-specific: for each SNP the lookup table says which
forward-strand nucleotide the manifest's allele A and allele B correspond to.
``AB -> FORWARD`` replaces the labels A/B with those nucleotides;
``FORWARD -> AB`` does the reverse. Missing calls stay missing.

Allele order in the output is alphabetical (``AG``, not ``GA``) so results
are deterministic; see :func:`soysnp_format.formats.format_genotype`.
"""

from . import formats
from .detect import detect_format, AB, FORWARD, INDETERMINATE, MIXED
from .formats import format_genotype
from .io import GenotypeData
from .lookup import load_lookup


class ConversionError(Exception):
    """Raised when conversion cannot proceed (unknown SNPs, bad alleles, ...)."""


def _allele_map(from_format, to_format, pair):
    a_fwd, b_fwd = pair
    if from_format == AB and to_format == FORWARD:
        return {"A": a_fwd, "B": b_fwd}
    if from_format == FORWARD and to_format == AB:
        if a_fwd == b_fwd:
            raise ConversionError(
                f"lookup row has identical forward alleles ({a_fwd}); cannot invert"
            )
        return {a_fwd: "A", b_fwd: "B"}
    raise ConversionError(f"unsupported conversion {from_format} -> {to_format}")


def convert_genotypes(data, from_format, to_format, lookup=None, chip="auto",
                      on_unknown_snp="error", on_bad_allele="error"):
    """Convert all genotypes in ``data`` from one encoding to another.

    ``from_format`` may be ``"AB"``, ``"FORWARD"`` or ``"auto"`` (detected).
    ``on_unknown_snp``: ``"error"`` (default) or ``"skip"`` SNPs absent from
    the lookup. ``on_bad_allele``: ``"error"`` (default) or ``"missing"``
    (genotypes with alleles outside the expected source alphabet become
    missing calls).
    Returns a new :class:`GenotypeData` in wide (SNP-rows) layout.
    """
    from_format = (from_format or "").upper()
    to_format = (to_format or "").upper()
    if from_format == "AUTO":
        result = detect_format(data, lookup=lookup, chip=chip)
        if result.verdict in (INDETERMINATE, MIXED):
            raise ConversionError(
                f"cannot auto-detect source format ({result.verdict}, "
                f"confidence {result.confidence}); pass --from explicitly."
            )
        from_format = result.verdict
    if from_format == to_format:
        raise ConversionError("source and target formats are identical; nothing to do")
    if from_format not in formats.FORMATS or to_format not in formats.FORMATS:
        raise ConversionError(
            f"unsupported conversion {from_format} -> {to_format}; "
            f"expected one of {formats.FORMATS}"
        )
    if lookup is None:
        lookup = load_lookup(chip)

    unknown_snps = []
    bad_alleles = []  # (snp, sample, genotype)
    out_calls = {}
    for snp in data.snps:
        if snp not in lookup:
            unknown_snps.append(snp)
            if on_unknown_snp == "skip":
                continue
            continue  # collected; raised below for "error"
        amap = _allele_map(from_format, to_format, lookup.alleles(snp))
        # Alleles each genotype may legitimately use: the A/B labels, or --
        # for forward input -- exactly the SNP's own forward allele pair.
        # (A forward-format genotype using any other nucleotides contradicts
        # the chip definition, e.g. a "CC" call on an A/G SNP.)
        allowed = {"A", "B"} if from_format == AB else set(lookup.alleles(snp))
        out_calls[snp] = {}
        for sample, gt in data.iter_genotypes(snp):
            if gt is None:
                out_calls[snp][sample] = None
                continue
            if gt[0] not in allowed or gt[1] not in allowed:
                bad_alleles.append((snp, sample, gt))
                if on_bad_allele == "missing":
                    out_calls[snp][sample] = None
                    continue
                continue  # collected; raised below for "error"
            out_calls[snp][sample] = (amap[gt[0]], amap[gt[1]])

    problems = []
    if unknown_snps and on_unknown_snp == "error":
        problems.append(
            f"{len(unknown_snps)} SNP(s) not in the {lookup.name} lookup table "
            f"(e.g. {unknown_snps[:5]}); use --lookup with a mapping for your chip "
            "or --on-unknown-snp skip"
        )
    if bad_alleles and on_bad_allele == "error":
        examples = [f"{snp}/{sample}={''.join(gt)}" for snp, sample, gt in bad_alleles[:5]]
        problems.append(
            f"{len(bad_alleles)} genotype(s) use alleles outside the {from_format} "
            f"alphabet (e.g. {examples}); use --on-bad-allele missing to blank them"
        )
    if problems:
        raise ConversionError("; ".join(problems))

    return GenotypeData(out_calls, data.samples, layout="wide",
                        source=f"{data.source} [{from_format}->{to_format}]")


def write_wide_matrix(data, path, delimiter="\t", missing="--"):
    """Write genotype data as a SNP-rows x sample-columns matrix."""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(delimiter.join(["SNP"] + data.samples) + "\n")
        for snp in data.snps:
            row = [snp]
            for sample in data.samples:
                gt = data.calls[snp].get(sample)
                row.append(format_genotype(*(gt if gt else (None, None)), missing=missing))
            fh.write(delimiter.join(row) + "\n")


def write_long_table(data, path, delimiter="\t", missing="--", encoding_name="Forward"):
    """Write genotype data in a Final-Report-like long table."""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(delimiter.join(
            ["SNP Name", "Sample ID",
             f"Allele1 - {encoding_name}", f"Allele2 - {encoding_name}"]) + "\n")
        for snp in data.snps:
            for sample in data.samples:
                gt = data.calls[snp].get(sample)
                if gt is None:
                    a1 = a2 = missing
                else:
                    a1, a2 = sorted(gt)
                fh.write(delimiter.join([snp, sample, a1, a2]) + "\n")
