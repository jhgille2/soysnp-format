# Building the bundled lookup tables

This directory holds the reproducible build for the SNP-specific Illumina A/B
→ forward-strand allele mappings shipped in `../src/soysnp_format/tables/`.

## Source files (downloaded 2026-09-28)

| File | Source | Size |
|---|---|---|
| `SoySNP50K_TableS1.xlsx` | Song et al. 2013, PLoS ONE, Table S1 — https://pmc.ncbi.nlm.nih.gov/articles/PMC3555945/ (supplementary file `pone.0054985.s001.xlsx`) | 10,707,943 bytes |
| `glyma.Wm82.gnm1.mrk.SoySNP50K.gff3.gz` | SoyBase data store — https://data.soybase.org/Glycine/max/markers/Wm82.gnm1.mrk.SoySNP50K/ | 722,189 bytes |
| `glyma.Wm82.gnm2.mrk.SoySNP50K.gff3.gz` | SoyBase data store — https://data.soybase.org/Glycine/max/markers/Wm82.gnm2.mrk.SoySNP50K/ | 800,327 bytes |
| `glyma.Wm82.gnm2.mrk.SoySNP6K.gff3.gz` | SoyBase data store — https://data.soybase.org/Glycine/max/markers/Wm82.gnm2.mrk.SoySNP6K/ | 81,266 bytes |

The `.xlsx` and the two gnm2 GFF3 files were fetched with a browser download
task; the gnm1 GFF3 was fetched with `curl` after the build showed Table S1
coordinates match the gnm1 assembly (see below). The input files are vendored
in this directory so the tables can be rebuilt without re-downloading.

## Method (`build_lookups.py`)

1. **Read Table S1** (60,800 SNPs): BARC SNP id, chromosome, coordinate,
   Williams 82 allele, alternative allele, and the 60 bp flanking sequence
   with the SNP shown as `[X/Y]`.
2. **Apply Illumina's published A/B allele rule** to the forward-strand
   alleles from the flanking sequence. Source: Illumina Technical Note
   *"TOP/BOT Strand and A/B Allele: A guide to Illumina's method for
   determining Strand and Allele for the GoldenGate and Infinium Assays"*
   (Pub. No. 370-2006-018, 26 Jun 2006):
   - SNP alleles containing **A** with partner C/G → strand TOP; Allele A = A.
   - SNP alleles containing **T** (no A) with partner C/G → Allele A = T
     (strand designated BOT; the A/B labels attach to the alleles on the
     submitted/forward strand — Table 1 of the note).
   - `[A/T]` and `[C/G]` SNPs are ambiguous and need "sequence walking" on the
     flanking sequence (implemented, with the decisive pair and walk distance
     recorded in the `derivation` column).
3. **Join to the gnm1 GFF3** on (chromosome, position) — Table S1 coordinates
   are on that assembly — to attach dbSNP ss IDs. All 60,800 join exactly and
   the GFF3 forward alleles agree with the Table S1 bracket alleles in every
   case, confirming the bracket alleles are forward-strand.
4. **Join gnm1 → gnm2 by ss ID** for current-assembly positions. Where gnm2
   forward alleles are the complement of the Table S1 alleles (assembly
   strand flip), the derived A/B labels are complemented; SNPs matching
   neither are excluded.
5. **BARCSoySNP6K**: the 5,989 SNPs in the 6K GFF3 are joined by ss ID to the
   50K-derived mapping and inherit it, cross-checked against the 6K GFF3
   alleles.

## Build results (2026-09-28 run)

- Table S1 records: 60,800. A/B derivation: **60,800 unambiguous, 0 sequence
  walks, 0 failures**. Allele-pair census: A/G 23,692; C/T 23,946; A/C 6,589;
  G/T 6,573; **A/T 0; C/G 0**. The chip contains no ambiguous SNP types, so
  every mapping follows directly from the allele pair — no flanking-sequence
  interpretation was needed for any SNP.
- gnm1 join: 60,800 exact position matches, 0 allele mismatches.
- gnm2 join: 60,556 allele-set matches, 0 complemented, 0 mismatched,
  244 ss IDs absent from the gnm2 file (these rows keep their gnm1 position,
  noted in the `notes` column).
- 6K: 5,989 inherited, 0 missing, 0 allele mismatches.
- Output: `../src/soysnp_format/tables/soysnp50k.csv` (60,800 rows),
  `../src/soysnp_format/tables/soysnp6k.csv` (5,989 rows).

## Caveats

- The mappings were **derived from the published Illumina rule**, not read from
  the official chip manifests (which are not publicly available). The rule is
  the same deterministic procedure Illumina uses to generate manifest A/B
  assignments, and the derivation was validated against the worked examples in
  the technical note itself (Tables 1–3).
- **Manifest validation (2026-09-28).** The BARCSoySNP6K mappings were checked
  SNP-by-SNP against the official `BARCSoySNP6k_11691901_A` manifest: for each
  of the 5,392 SNPs present in both, the manifest's `SNP` column `[X/Y]`
  (allele A / allele B on the assay's `IlmnStrand`) equals the table's
  `(a_forward, b_forward)` or its strand complement, with the A/B labels never
  swapped. Result: **5,392 agree, 0 disagree** (11 manifest SNPs have no
  counterpart in this table; 597 table SNPs are not on that manifest version).
  Re-run with `python validate_against_manifest.py MANIFEST.csv --table
  ../src/soysnp_format/tables/soysnp6k.csv`. The SoySNP50K table uses the same
  derivation procedure but has not been checked against its manifest.
- Positions are Wm82.gnm2 (falling back to gnm1 for the 244 unmapped SNPs);
  the A/B → forward *allele* mapping does not depend on the assembly.
