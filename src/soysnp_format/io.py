"""Read GenomeStudio exports and generic SNP array matrices.

Two input layouts are supported:

Final Report (``final-report``)
    GenomeStudio's "Final Report" text export. It starts with ``[Header]`` and
    ``[Data]`` sections; the data section is a long table with columns
    ``SNP Name``, ``Sample ID``, ``Allele1 - <enc>``, ``Allele2 - <enc>``
    where ``<enc>`` is the export encoding GenomeStudio used (e.g. ``AB``,
    ``Forward``, ``Top``). The encoding named in the column headers is kept
    as a *declared-format hint* -- it is checked against the data, not trusted
    blindly.

Wide matrix (``wide``)
    A delimited table with SNP identifiers on one axis and sample identifiers
    on the other. Orientation (SNPs-as-rows vs SNPs-as-columns) is inferred
    from the bundled chip SNP lists when possible, otherwise from simple
    heuristics (see :func:`read_genotypes`).

The parsed result is a :class:`GenotypeData`: ``snp -> {sample: genotype}``
with genotypes as 2-tuples from :func:`soysnp_format.formats.normalize_genotype`
(``None`` for missing calls).
"""

import csv
import re

from . import formats
from .formats import normalize_genotype

#: Maps GenomeStudio "Allele1 - <enc>" header suffixes to canonical formats.
DECLARED_FORMAT_HINTS = {
    "AB": formats.AB,
    "A/B": formats.AB,
    "FORWARD": formats.FORWARD,
    "PLUS": formats.FORWARD,   # plus strand == forward strand for these exports
    "TOP": None,               # recognized, but not translated by this tool
    "BOTTOM": None,
    "BOT": None,
    "MINUS": None,
    "DESIGN": None,
}


class GenotypeData:
    """Parsed genotype matrix: snp_id -> {sample_id: genotype tuple or None}."""

    def __init__(self, calls, samples, declared_format=None, layout=None, source=None):
        self.calls = calls            # dict snp -> dict sample -> (a1, a2) | None
        self.samples = list(samples)  # sample ids in file order
        self.snps = list(calls)       # snp ids in file order
        self.declared_format = declared_format  # format hint from file headers, if any
        self.layout = layout          # 'final-report' or 'wide'
        self.source = source

    def __len__(self):
        return len(self.snps)

    def iter_genotypes(self, snp):
        """Yield (sample, genotype) pairs for one SNP."""
        return self.calls[snp].items()


def _sniff_delimiter(sample_text):
    try:
        dialect = csv.Sniffer().sniff(sample_text, delimiters=",\t;|")
        return dialect.delimiter
    except csv.Error:
        # fall back to whichever of tab/comma appears more often
        return "\t" if sample_text.count("\t") >= sample_text.count(",") else ","


def _declared_format_from_header(columns):
    """Extract GenomeStudio allele-encoding hint from 'Allele1 - <enc>' headers."""
    for col in columns:
        m = re.match(r"(?i)^allele\s*[12]\s*-\s*(.+?)\s*$", col.strip())
        if m:
            enc = m.group(1).strip().upper()
            return DECLARED_FORMAT_HINTS.get(enc, "UNKNOWN:" + enc)
    return None


def _try_genotype(token):
    """Like :func:`normalize_genotype`, but returns ``None`` for garbage."""
    try:
        return normalize_genotype(token)
    except ValueError:
        return None


def _read_final_report_matrix(path, lines, data_start, delim, columns):
    """Parse a matrix-style Final Report ``[Data]`` section.

    GenomeStudio can export the ``[Data]`` section as a SNP x sample matrix
    instead of the long SNP Name / Sample ID / Allele1 / Allele2 layout: the
    first data row holds sample IDs (with an empty stub cell), and each
    following row holds one SNP ID plus one genotype call per sample
    (e.g. ``AA``/``AB``/``BB``).

    Returns a :class:`GenotypeData` or ``None`` if the header row does not
    look like a sample-ID row.
    """
    if len(columns) < 2:
        return None
    stub, sample_cells = columns[0].strip(), columns[1:]
    # the stub cell must be empty or a stub label, not a SNP id or a call
    if stub and stub.lower() not in ("snp", "snp name", "snpname", "marker",
                                     "locus", "sample", "sample id"):
        return None
    samples = [c.strip() for c in sample_cells if c.strip()]
    if not samples:
        return None
    # guard: if the "sample" cells are mostly genotype calls, this row is a
    # data row, not a header row -> not a matrix-style report
    geno_like = sum(1 for c in sample_cells
                    if _try_genotype(c.strip()) is not None)
    if geno_like > len(sample_cells) / 2:
        return None

    calls = {}
    reader = csv.reader(lines[data_start + 1:], delimiter=delim)
    for row in reader:
        if not row or not any(c.strip() for c in row):
            continue
        snp = row[0].strip()
        if not snp or snp.startswith("["):
            continue
        calls[snp] = {}
        for sample, token in zip(samples, row[1:]):
            calls[snp][sample] = normalize_genotype(token.strip())
    if not calls:
        return None
    return GenotypeData(calls, samples, declared_format=None,
                        layout="final-report", source=str(path))


def read_final_report(path, encoding="utf-8-sig"):
    """Parse a GenomeStudio Final Report file.

    Handles both the long SNP Name / Sample ID / Allele1 / Allele2 layout
    and the matrix-style ``[Data]`` section (sample IDs across the first
    row, one SNP per row).
    """
    with open(path, "r", encoding=encoding, newline="") as fh:
        text = fh.read()
    lines = text.splitlines()
    # find the [Data] section
    data_start = None
    for i, line in enumerate(lines):
        if line.strip().lower() == "[data]":
            data_start = i + 1
            break
    if data_start is None:
        raise ValueError(f"{path}: not a GenomeStudio Final Report ([Data] section not found)")
    header = lines[data_start].rstrip("\n")
    delim = _sniff_delimiter(header + "\n" + "\n".join(lines[data_start + 1:data_start + 6]))
    columns = next(csv.reader([header], delimiter=delim))
    col_index = {c.strip().lower(): i for i, c in enumerate(columns)}

    def find(*names):
        for n in names:
            if n in col_index:
                return col_index[n]
        return None

    snp_col = find("snp name", "snpname", "snp", "locus")
    sample_col = find("sample id", "sampleid", "sample", "sample name", "samplename")
    a1_col = find("allele1 - ab", "allele1-ab")
    a2_col = find("allele2 - ab", "allele2-ab")
    if a1_col is None:  # fall back to any Allele1/Allele2 pair
        a1_col = next((i for i, c in enumerate(columns)
                       if re.match(r"(?i)^allele\s*1\b", c.strip())), None)
        a2_col = next((i for i, c in enumerate(columns)
                       if re.match(r"(?i)^allele\s*2\b", c.strip())), None)
    if snp_col is None or sample_col is None or a1_col is None or a2_col is None:
        matrix = _read_final_report_matrix(
            path, lines, data_start, delim, columns)
        if matrix is not None:
            return matrix
        raise ValueError(
            f"{path}: Final Report data header lacks SNP Name / Sample ID / "
            f"Allele1 / Allele2 columns (saw: {columns[:8]}...)"
        )
    declared = _declared_format_from_header(columns)

    calls = {}
    samples = []
    seen_samples = set()
    reader = csv.reader(lines[data_start + 1:], delimiter=delim)
    for row in reader:
        if not row or not any(c.strip() for c in row):
            continue
        if len(row) <= max(snp_col, sample_col, a1_col, a2_col):
            continue
        snp = row[snp_col].strip()
        sample = row[sample_col].strip()
        if not snp or not sample:
            continue
        gt = normalize_genotype(row[a1_col].strip() + row[a2_col].strip())
        if sample not in seen_samples:
            seen_samples.add(sample)
            samples.append(sample)
        calls.setdefault(snp, {})[sample] = gt
    return GenotypeData(calls, samples, declared_format=declared,
                        layout="final-report", source=str(path))


def _looks_like_snp_id(token, snp_sets):
    token = token.strip()
    return any(token in s for s in snp_sets if s)


def read_wide_matrix(path, snp_sets=(), encoding="utf-8-sig"):
    """Parse a wide SNP x sample (or sample x SNP) matrix.

    ``snp_sets`` is an optional iterable of sets of known SNP ids used to
    decide the matrix orientation. Without it, orientation is guessed: the
    axis whose labels look like SNP ids (``ss...``/``ss715...`` style or
    simply the longer axis labels) is treated as SNPs.
    """
    with open(path, "r", encoding=encoding, newline="") as fh:
        sample_text = fh.read(65536)
        fh.seek(0)
        delim = _sniff_delimiter(sample_text)
        reader = csv.reader(fh, delimiter=delim)
        rows = [[c.strip() for c in row] for row in reader if any(c.strip() for c in row)]
    # tolerate a GenomeStudio [Header]/[Data] preamble: parse from [Data]
    if rows and rows[0] and rows[0][0].startswith("["):
        for i, r in enumerate(rows):
            if r and r[0].strip().lower() == "[data]":
                rows = rows[i + 1:]
                break
    if len(rows) < 2 or len(rows[0]) < 2:
        raise ValueError(f"{path}: wide matrix needs at least 2 rows and 2 columns")

    header = rows[0]
    first_col = [r[0] for r in rows[1:]]
    header_ids = [h for h in header[1:] if h]
    first_col_ids = [c for c in first_col if c]

    snps_are_rows = None
    if snp_sets:
        header_hits = sum(1 for h in header_ids if _looks_like_snp_id(h, snp_sets))
        col_hits = sum(1 for c in first_col_ids if _looks_like_snp_id(c, snp_sets))
        if header_hits or col_hits:
            # SNPs on the axis matching more known chip SNP ids
            snps_are_rows = col_hits >= header_hits

    if snps_are_rows is None:
        # heuristic: SNP ids often look like ss715... ; otherwise assume the
        # longer axis holds the SNPs (arrays type many more SNPs than samples)
        def snpish(tok):
            return bool(re.match(r"(?i)^(ss|rs|tm|sb)\d+$", tok))
        header_hits = sum(1 for h in header_ids if snpish(h))
        col_hits = sum(1 for c in first_col_ids if snpish(c))
        if header_hits or col_hits:
            snps_are_rows = col_hits >= header_hits
        else:
            snps_are_rows = len(first_col_ids) >= len(header_ids)

    calls = {}
    if snps_are_rows:
        samples = header_ids
        for row in rows[1:]:
            snp = row[0]
            if not snp:
                continue
            calls[snp] = {}
            for sample, token in zip(samples, row[1:]):
                calls[snp][sample] = normalize_genotype(token)
    else:
        samples = first_col_ids
        snp_ids = header_ids
        for snp in snp_ids:
            calls[snp] = {}
        for row in rows[1:]:
            sample = row[0]
            if not sample:
                continue
            for snp, token in zip(snp_ids, row[1:]):
                calls[snp][sample] = normalize_genotype(token)
    return GenotypeData(calls, samples, layout="wide", source=str(path))


def read_genotypes(path, layout="auto", snp_sets=(), encoding="utf-8-sig"):
    """Read ``path`` and return a :class:`GenotypeData`.

    ``layout`` is ``auto`` (default), ``final-report`` or ``wide``.
    ``snp_sets`` optionally supplies known SNP-id sets to orient wide matrices.
    """
    with open(path, "r", encoding=encoding) as fh:
        head = fh.read(4096)
    head_low = head.lower()
    is_final_report = "[data]" in head_low and "[header]" in head_low
    if layout == "auto":
        layout = "final-report" if is_final_report else "wide"
    if layout == "final-report":
        return read_final_report(path, encoding=encoding)
    if layout == "wide":
        return read_wide_matrix(path, snp_sets=snp_sets, encoding=encoding)
    raise ValueError(f"unknown layout {layout!r}: expected 'auto', 'final-report' or 'wide'")
