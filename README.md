# soysnp-format

Detect and translate allele-encoding formats (Illumina A/B vs forward strand)
in soybean SNP array exports.

Legacy soybean genotype files are ambiguous: a file full of `AA`/`AB`/`BB`
calls uses Illumina's A/B encoding, while one full of `AG`/`CT` calls uses
forward-strand nucleotides — but a file containing only `AA` calls could be
either. `soysnp-format` inspects the minimum information needed to tell them
apart, reports uncertainty honestly, and converts between encodings using
SNP-specific lookup tables.

Bundled chip tables: **SoySNP50K** (60,800 SNPs) and **BARCSoySNP6K**
(5,989 SNPs), with both BARC SNP ids and dbSNP ss ids recognized. See
[`build/README.md`](build/README.md) for how the tables were derived and
their provenance.

## Installation

```bash
pip install soysnp-format
```

or from source:

```bash
git clone https://github.com/<you>/soysnp-format.git
cd soysnp-format
pip install -e .
```

Requires Python ≥ 3.9. No other dependencies.

## Usage

```bash
# What encoding is this file in?
soysnp-format detect genotypes.tsv

# Convert A/B calls to forward-strand nucleotides
soysnp-format convert genotypes.tsv --from AB --to FORWARD -o genotypes_forward.tsv

# Convert with auto-detection of the source encoding
soysnp-format convert genotypes.tsv --from auto --to AB -o genotypes_ab.tsv

# Which chip tables are bundled?
soysnp-format chips
```

Accepted inputs: GenomeStudio Final Report exports (`[Header]`/`[Data]`
sections) and generic wide matrices (SNPs × samples, either orientation).
Use `--layout` to force one, `--chip` to restrict to one bundled table, or
`--lookup mychip.csv` for a custom chip table (`snp,a_forward,b_forward`).

### Detection output

```
Likely format: AB (confidence: moderate)
SNPs examined: 3
  with A/B-only evidence (a 'B' call seen): 3
  with forward-only evidence (C/G/T seen):  0
  ambiguous (only 'A' or no calls):         0
  conflicting (both B and C/G/T seen):      0
SNP ids matched against soysnp50k+soysnp6k: 3/3
```

Possible verdicts: `AB`, `FORWARD`, `MIXED` (different SNPs show different
encodings — often a sign of a malformed file), and `INDETERMINATE` (no
decisive evidence, e.g. only `AA` calls). Confidence is `high`, `moderate`,
or `low` depending on how much decisive evidence was seen. `--json` emits
machine-readable output.

### Conversion

Conversion is SNP-specific: each SNP's A/B alleles are translated through
that SNP's own forward-strand allele pair from the lookup table. SNPs missing
from the table abort with an error by default (`--on-unknown-snp skip` to
drop them instead); genotypes outside the source alphabet do the same
(`--on-bad-allele missing` to emit missing calls instead). Output layouts:
`wide` (SNP-rows matrix) or `long` (Final-Report-like table). Output alleles
are written in alphabetical order (`AG`, not `GA`).

## Web app

A static browser version lives in `docs/` (served via GitHub Pages): drop a
genotype file onto the page to detect its encoding and convert it, with no
installation and no data leaving your computer. The page embeds the same
detection/conversion logic as the Python package (`docs/soysnp.js`, tested
headlessly in `node`) and loads the bundled SNP lookup table (`docs/assets/lookups.json`)
on first use.

## As a library

```python
from soysnp_format import read_genotypes, detect_format, convert_genotypes
from soysnp_format.lookup import load_lookup

lookup = load_lookup("soysnp50k")          # or "soysnp6k", or "auto"
data = read_genotypes("genotypes.tsv")
print(detect_format(data, lookup=lookup).summary())
fwd = convert_genotypes(data, "AB", "FORWARD", lookup=lookup)
```

## What "forward" means here

The bundled tables map Illumina A/B alleles to forward-strand nucleotides as
Illumina's GenomeStudio "Forward" encoding defines them — i.e. the allele
labels applied to the sequence as submitted for the assay design. This is
verified to agree with the forward strand of the Williams 82 assemblies for
all 60,800 SoySNP50K SNPs (see `build/README.md`), but it is not a claim
about any other assembly or strand convention.

## Deriving the bundled tables

`build/build_lookups.py` regenerates `src/soysnp_format/tables/` from the
downloaded sources in `build/` (Table S1 + SoyBase GFF3 marker files).
Summary of results and caveats: [`build/README.md`](build/README.md).

## License

MIT. See [LICENSE](LICENSE).
