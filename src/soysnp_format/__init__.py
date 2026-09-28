"""soysnp_format: detect and translate allele-encoding formats in soybean SNP array exports.

Supported formats:
    AB       Illumina A/B encoding as exported by GenomeStudio (genotypes AA/AB/BB,
             where A and B are the manifest design alleles).
    FORWARD  Forward-strand nucleotide encoding (genotypes such as AA/AG/GG), the
             encoding used by most public repositories (e.g. SoyBase).

Bundled chip lookup tables cover the SoySNP50K iSelect BeadChip and the
BARCSoySNP6K assay (a ~6K subset of SoySNP50K SNPs, so both chips share one
allele-mapping table keyed by SNP name).
"""

from .detect import detect_format, DetectionResult
from .convert import convert_genotypes
from .io import read_genotypes, GenotypeData
from .lookup import ChipLookup, load_lookup

__all__ = [
    "detect_format",
    "DetectionResult",
    "convert_genotypes",
    "read_genotypes",
    "GenotypeData",
    "ChipLookup",
    "load_lookup",
]

__version__ = "0.1.0"
