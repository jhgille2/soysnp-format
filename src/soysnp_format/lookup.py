"""Per-SNP allele lookup tables: the bridge between A/B and FORWARD encodings.

A lookup table maps each SNP name to the forward-strand nucleotides that
correspond to the manifest's allele A and allele B::

    snp,a_forward,b_forward

For example, a row ``BARC_1.01_Gm01_2033_G_A,A,G`` means: in A/B encoding this
SNP's genotypes are AA/AB/BB, and on the forward strand allele A is ``A`` and
allele B is ``G`` (so AB-format ``AB`` == forward-format ``AG``).

Bundled tables
--------------
``tables/soysnp50k.csv`` and ``tables/soysnp6k.csv`` cover the SoySNP50K
iSelect BeadChip (60,800 SNPs) and the BARCSoySNP6K assay (5,989 SNPs). Every
BARCSoySNP6K SNP is also a SoySNP50K SNP, so the two tables share the same
allele mappings; the 6K table is the subset of rows for 6K SNPs.

The mappings were derived from the published Illumina A/B allele rule applied
to the SoySNP50K Table S1 sequences (Song et al. 2013), joined to SoyBase
Wm82 GFF3 marker files for ss IDs and current positions. Rows are keyed by
the BARC SNP id (e.g. ``BARC_1.01_Gm01_2033_G_A``) with the dbSNP ss ID
(e.g. ``ss715578672``) and the manifest-style short name (e.g.
``Gm01_2033_G_A``, as used in GenomeStudio exports) as aliases, so genotype
files using any of these naming schemes match. See ``build/README.md`` for
the full derivation, provenance, and validation status.

Custom tables
-------------
For any other chip, build a CSV with the columns above (e.g. from the chip
manifest: for each SNP record which forward-strand base the manifest's
allele A and allele B correspond to) and pass it with ``--lookup``.
An optional ``ss_id`` column is indexed as an alias.
"""

import csv
import os

TABLE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tables")

BUNDLED_CHIPS = ("soysnp50k", "soysnp6k")

#: Header names accepted for the SNP identifier column of a lookup CSV.
SNP_COLUMNS = ("snp", "snp_name", "snpname", "marker", "marker_name", "name", "ilmnid")


class ChipLookup:
    """Forward-strand allele mapping for one chip: snp -> (a_forward, b_forward).

    ``aliases`` maps alternative SNP ids (e.g. ss IDs) to primary SNP ids;
    membership tests, :meth:`alleles` and :meth:`snps` resolve them.
    """

    def __init__(self, mapping, name="custom", aliases=None):
        # mapping: dict snp -> (a_forward, b_forward), alleles uppercase single letters
        self.mapping = dict(mapping)
        self.name = name
        self.aliases = dict(aliases) if aliases else {}

    def _primary(self, snp):
        return self.aliases.get(snp, snp)

    def __contains__(self, snp):
        return self._primary(snp) in self.mapping

    def __len__(self):
        return len(self.mapping)

    def snps(self):
        return set(self.mapping) | set(self.aliases)

    def alleles(self, snp):
        """Return (a_forward, b_forward) for ``snp`` (or an alias); raises KeyError if absent."""
        return self.mapping[self._primary(snp)]

    def forward_pair(self, snp):
        """Return the unordered forward allele set {a_forward, b_forward}."""
        a, b = self.alleles(snp)
        return {a, b}


def _find_column(columns, candidates):
    lowered = {c.strip().lower(): c for c in columns}
    for cand in candidates:
        if cand in lowered:
            return lowered[cand]
    return None


def read_lookup_csv(path, name=None):
    """Read a lookup CSV with columns ``snp,a_forward,b_forward`` (see module docstring).

    An ``ss_id`` column, when present, is indexed as an alias for the primary
    SNP id. Extra columns are ignored.
    """
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: empty lookup file")
        snp_col = _find_column(reader.fieldnames, SNP_COLUMNS)
        a_col = _find_column(reader.fieldnames,
                             ("a_forward", "allelea_forward", "allele_a_forward", "forward_a"))
        b_col = _find_column(reader.fieldnames,
                             ("b_forward", "alleleb_forward", "allele_b_forward", "forward_b"))
        ss_col = _find_column(reader.fieldnames, ("ss_id", "ssid", "dbsnp"))
        if snp_col is None or a_col is None or b_col is None:
            raise ValueError(
                f"{path}: lookup CSV needs snp, a_forward and b_forward columns "
                f"(saw: {reader.fieldnames})"
            )
        mapping = {}
        aliases = {}
        for row in reader:
            snp = (row.get(snp_col) or "").strip()
            a = (row.get(a_col) or "").strip().upper()
            b = (row.get(b_col) or "").strip().upper()
            if not snp or len(a) != 1 or len(b) != 1:
                continue
            if a not in "ACGT" or b not in "ACGT":
                raise ValueError(f"{path}: SNP {snp}: forward alleles must be A/C/G/T")
            mapping[snp] = (a, b)
            # alias: manifest/GenomeStudio exports use the BARC id without the
            # 'BARC_1.01_' prefix (e.g. 'Gm01_1013695_A_G')
            if snp.startswith("BARC_1.01_"):
                short = snp[len("BARC_1.01_"):]
                if short and short not in mapping and short not in aliases:
                    aliases[short] = snp
            if ss_col:
                alias = (row.get(ss_col) or "").strip()
                if alias and alias not in mapping and alias not in aliases:
                    aliases[alias] = snp
    return ChipLookup(mapping, name=name or os.path.basename(path), aliases=aliases)


def load_lookup(chip="auto", path=None):
    """Load a :class:`ChipLookup`.

    ``chip``: ``"soysnp50k"``, ``"soysnp6k"`` or ``"auto"`` (both bundled
    tables merged; 6K SNPs are a subset of 50K so mappings agree).
    ``path``: optional custom lookup CSV, used instead of the bundled tables.
    """
    if path:
        return read_lookup_csv(path, name="custom")
    chip = (chip or "auto").lower()
    if chip == "auto":
        merged: dict = {}
        merged_aliases: dict = {}
        for c in BUNDLED_CHIPS:
            table = load_lookup(c)
            merged.update(table.mapping)
            merged_aliases.update(table.aliases)
        return ChipLookup(merged, name="soysnp50k+soysnp6k", aliases=merged_aliases)
    if chip not in BUNDLED_CHIPS:
        raise ValueError(f"unknown chip {chip!r}: expected one of {BUNDLED_CHIPS} or a --lookup file")
    table_path = os.path.join(TABLE_DIR, f"{chip}.csv")
    if not os.path.exists(table_path):
        raise FileNotFoundError(
            f"bundled lookup table not found: {table_path}. "
            "See build/README.md for how the tables are generated."
        )
    return read_lookup_csv(table_path, name=chip)


def bundled_snp_sets():
    """Return {chip: set_of_snp_ids} for the bundled tables (used for orientation)."""
    sets = {}
    for chip in BUNDLED_CHIPS:
        try:
            sets[chip] = load_lookup(chip).snps()
        except FileNotFoundError:
            sets[chip] = set()
    return sets
