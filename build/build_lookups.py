#!/usr/bin/env python3
"""Derive SNP-specific Illumina A/B -> forward-strand mappings for SoySNP50K / BARCSoySNP6K.

Method (documented in ../BUILD.md):
  1. Read SoySNP50K Table S1 (Song et al. 2013): per-SNP BARC id, chromosome,
     position, Wm82 allele, alternative allele, and the 60 bp flanking sequence
     with the SNP shown as [X/Y].
  2. Apply Illumina's published A/B allele rule (Illumina Technical Note,
     "TOP/BOT Strand and A/B Allele: A guide to Illumina's method for
     determining Strand and Allele for the GoldenGate and Infinium Assays",
     Pub. No. 370-2006-018):
       - SNP alleles containing A with partner C/G  -> TOP; Allele A = A.
       - SNP alleles containing T (no A) with partner C/G -> Allele A = T
         (strand designated BOT, but the A/B labels attach to the alleles on
         the submitted/forward strand, Table 1 of the note).
       - [A/T] and [C/G] SNPs are ambiguous: "sequence walking" outward from
         the SNP finds the first unambiguous pair (A/G, A/C, T/C, T/G) at
         symmetric positions n-k/n+k. If the A or T of that pair is 5' of the
         SNP the sequence is TOP, else BOT. Then:
             [A/T] TOP -> A=A,B=T ; BOT -> A=T,B=A
             [C/G] TOP -> A=C,B=G ; BOT -> A=G,B=C
     The rule is strand-consistent (Table 3 of the note): applying it to the
     forward-strand flanking sequence labels the same physical alleles as
     applying it to the reverse complement.
  3. Join to the SoyBase Wm82.gnm1 GFF3 marker file on (chromosome, position)
     -- the assembly Table S1 coordinates match -- to attach ss IDs, then join
     gnm1 -> gnm2 by ss ID for current positions. Where gnm2 forward alleles
     are the complement of the Table S1 alleles (strand flip between
     assemblies), the derived A/B labels are complemented accordingly; SNPs
     whose allele sets match neither are flagged and excluded.
  4. BARCSoySNP6K SNPs are the subset of SoySNP50K SNPs present in the 6K GFF3
     (joined by ss ID); they inherit the 50K-derived mapping, cross-checked
     against the 6K GFF3 alleles.

Outputs (written to ../src/soysnp_format/tables/):
  - soysnp50k.csv
  - soysnp6k.csv

This script is committed for reproducibility; re-running it requires the three
downloaded source files in this directory (see ../BUILD.md for URLs).
"""

from __future__ import annotations

import csv
import gzip
import re
import sys
from pathlib import Path

BUILD_DIR = Path(__file__).resolve().parent
TABLE_DIR = BUILD_DIR.parent / "src" / "soysnp_format" / "tables"

TABLE_S1 = BUILD_DIR / "SoySNP50K_TableS1.xlsx"
GFF1K = BUILD_DIR / "glyma.Wm82.gnm1.mrk.SoySNP50K.gff3.gz"
GFF50K = BUILD_DIR / "glyma.Wm82.gnm2.mrk.SoySNP50K.gff3.gz"
GFF6K = BUILD_DIR / "glyma.Wm82.gnm2.mrk.SoySNP6K.gff3.gz"

FLANK_RE = re.compile(r"^([ACGTNacgtn]+)\[([ACGT])/([ACGT])\]([ACGTNacgtn]+)$")
COMP = str.maketrans("ACGTN", "TGCAN")


def complement(base: str) -> str:
    return base.translate(COMP)


# ---------------------------------------------------------------------------
# Illumina A/B rule
# ---------------------------------------------------------------------------

UNAMBIGUOUS_PAIRS = {
    frozenset("AG"), frozenset("AC"), frozenset("TC"), frozenset("TG"),
}


def sequence_walk(left: str, right: str) -> tuple[str, str] | None:
    """Return (strand, detail) via Illumina sequence walking, or None.

    ``left``/``right`` are the 5'/3' flanking sequences (5'->3' orientation).
    Walks k = 1..len outward; the first unambiguous pair (A/G, A/C, T/C, T/G)
    at (n-k, n+k) decides: A-or-T on the 5' side -> TOP, on the 3' side -> BOT.
    """
    n = min(len(left), len(right))
    for k in range(1, n + 1):
        l, r = left[-k].upper(), right[k - 1].upper()
        pair = frozenset((l, r))
        if pair in UNAMBIGUOUS_PAIRS and l != r:
            at = "A" if "A" in pair else "T"  # the A-or-T member of the pair
            side = "5p" if (l == at) else "3p"
            strand = "TOP" if side == "5p" else "BOT"
            return strand, f"walk:k={k}:pair={l}/{r}@{side}->{strand}"
    return None


def illumina_ab(snp_x: str, snp_y: str, left: str, right: str):
    """Return (a_forward, b_forward, derivation) for forward alleles X/Y.

    Raises ValueError if the SNP is not a biallelic ACGT SNP or the walk fails.
    """
    alleles = frozenset((snp_x, snp_y))
    if alleles == frozenset("AT"):
        walked = sequence_walk(left, right)
        if walked is None:
            raise ValueError("sequence walk found no unambiguous pair")
        strand, detail = walked
        if strand == "TOP":
            return "A", "T", detail + ":A/T:TOP"
        return "T", "A", detail + ":A/T:BOT"
    if alleles == frozenset("CG"):
        walked = sequence_walk(left, right)
        if walked is None:
            raise ValueError("sequence walk found no unambiguous pair")
        strand, detail = walked
        if strand == "TOP":
            return "C", "G", detail + ":C/G:TOP"
        return "G", "C", detail + ":C/G:BOT"
    if "A" in alleles:  # A/G or A/C -> TOP, Allele A = A
        other = next(iter(alleles - {"A"}))
        return "A", other, f"unambiguous:A/{other}:TOP"
    if "T" in alleles:  # T/C or T/G -> Allele A = T (Table 1: BOT, labels on strand)
        other = next(iter(alleles - {"T"}))
        return "T", other, f"unambiguous:T/{other}:BOT-design"
    raise ValueError(f"unexpected allele set {sorted(alleles)}")


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def read_table_s1(path: Path):
    import openpyxl

    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[wb.sheetnames[0]]
    rows = list(ws.iter_rows(values_only=True))
    header = [str(c).strip() for c in rows[1]]
    idx = {name: i for i, name in enumerate(header)}
    records = []
    for row in rows[2:]:
        if not row[idx["SNP ID"]]:
            continue
        records.append({
            "barcode_id": str(row[idx["SNP ID"]]).strip(),
            "chromosome": str(row[idx["Chromosome"]]).strip(),
            "position": int(row[idx["Coordinate"]]),
            "w82_allele": str(row[idx["Williams 82 allele"]]).strip().upper(),
            "alt_allele": str(row[idx["Alternative allele"]]).strip().upper(),
            "flank": str(row[idx["60bp sequence flanking SNP"]]).strip(),
        })
    return records


def read_gff3(path: Path):
    """Return ({(chrom, pos): record}, {ss_id: record})."""
    by_pos: dict[tuple[str, int], dict] = {}
    by_ss: dict[str, dict] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos = fields[0], int(fields[3])
            attrs = dict(kv.split("=", 1) for kv in fields[8].split(";") if "=" in kv)
            ss_id = attrs.get("Name", "")
            alleles = tuple(a.strip().upper() for a in attrs.get("alleles", "").split("/"))
            rec = {"chrom": chrom, "pos": pos, "ss_id": ss_id, "alleles": alleles}
            by_pos[(chrom, pos)] = rec
            if ss_id:
                by_ss[ss_id] = rec
    return by_pos, by_ss


def short_chrom(chrom: str) -> str:
    """'glyma.Wm82.gnm2.Gm01' -> 'Gm01'."""
    return chrom.split(".")[-1]


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

OUT_FIELDS = ["snp", "ss_id", "chromosome", "position",
              "a_forward", "b_forward", "derivation", "alleles_gnm2", "notes"]


def main() -> int:
    for p in (TABLE_S1, GFF1K, GFF50K, GFF6K):
        if not p.exists():
            print(f"missing input: {p}", file=sys.stderr)
            return 2

    print("reading Table S1 ...", flush=True)
    s1 = read_table_s1(TABLE_S1)
    print(f"  {len(s1)} SNP records", flush=True)

    print("reading GFF3 files ...", flush=True)
    gff1_pos, gff1_ss = read_gff3(GFF1K)
    gff50_pos, gff50_ss = read_gff3(GFF50K)
    gff6_pos, gff6_ss = read_gff3(GFF6K)
    print(f"  gnm1 50K GFF3: {len(gff1_pos)} markers; "
          f"gnm2 50K GFF3: {len(gff50_pos)} markers; "
          f"gnm2 6K GFF3: {len(gff6_pos)} markers", flush=True)

    def find_gff1(chrom: str, pos: int):
        for cand in ((chrom, pos),
                     ("glyma.Wm82.gnm1." + chrom, pos),
                     ("glyma.Wm82.gnm2." + chrom, pos)):
            for table in (gff1_pos, gff50_pos):
                if cand in table:
                    return table[cand]
        return None

    derived = []       # per-SNP derived records
    stats = {"unambiguous": 0, "walk": 0, "walk_fail": 0, "non_acgt": 0,
             "gff1_match": 0, "gff1_allele_mismatch": 0, "gff1_absent": 0,
             "gnm2_match": 0, "gnm2_complement": 0, "gnm2_mismatch": 0,
             "gnm2_absent": 0}
    walk_ks: list[int] = []

    for rec in s1:
        m = FLANK_RE.match(rec["flank"])
        if not m:
            stats["non_acgt"] += 1
            continue
        left, x, y, right = m.group(1).upper(), m.group(2), m.group(3), m.group(4).upper()
        try:
            a_fwd, b_fwd, derivation = illumina_ab(x, y, left, right)
        except ValueError:
            stats["walk_fail"] += 1
            continue
        if derivation.startswith("unambiguous"):
            stats["unambiguous"] += 1
        else:
            stats["walk"] += 1
            walk_ks.append(int(re.search(r"k=(\d+)", derivation).group(1)))

        # join to gnm1 GFF3 on (chromosome, position): same assembly as Table S1
        notes = []
        g1 = find_gff1(rec["chromosome"], rec["position"])
        if g1 is None:
            stats["gff1_absent"] += 1
            notes.append("no gnm1 GFF3 position match; EXCLUDED")
            continue
        s1_set = frozenset((x, y))
        if frozenset(g1["alleles"]) != s1_set:
            stats["gff1_allele_mismatch"] += 1
            notes.append(f"gnm1 alleles {g1['alleles']} != Table S1 {x}/{y}; EXCLUDED")
            continue
        stats["gff1_match"] += 1
        ss_id = g1["ss_id"]

        # gnm1 -> gnm2 by ss ID for current-assembly position/alleles
        gnm2_pos, gnm2_alleles = "", ""
        g2 = gff50_ss.get(ss_id)
        if g2 is None:
            stats["gnm2_absent"] += 1
            notes.append("ss ID absent from gnm2 GFF3; position shown is Table S1/gnm1")
            gnm2_pos = str(g1["pos"])
        else:
            gff_set = frozenset(g2["alleles"])
            if gff_set == s1_set:
                stats["gnm2_match"] += 1
            elif gff_set == frozenset(complement(b) for b in s1_set):
                stats["gnm2_complement"] += 1
                a_fwd, b_fwd = complement(a_fwd), complement(b_fwd)
                notes.append("gnm2 forward strand complemented vs Table S1; A/B complemented")
            else:
                stats["gnm2_mismatch"] += 1
                notes.append(f"gnm2 alleles {g2['alleles']} match neither {x}/{y} nor complement; EXCLUDED")
                continue
            gnm2_pos = str(g2["pos"])
            gnm2_alleles = "/".join(g2["alleles"])

        derived.append({
            "snp": rec["barcode_id"],
            "ss_id": ss_id,
            "chromosome": rec["chromosome"],
            "position": gnm2_pos,
            "a_forward": a_fwd,
            "b_forward": b_fwd,
            "derivation": derivation,
            "alleles_gnm2": gnm2_alleles,
            "notes": "; ".join(notes),
        })

    # 6K: subset by ss_id, cross-check alleles
    k6, k6_missing, k6_allele_mismatch = [], 0, 0
    by_ss_derived = {d["ss_id"]: d for d in derived if d["ss_id"]}
    for ss_id, grec in gff6_ss.items():
        d = by_ss_derived.get(ss_id)
        if d is None:
            k6_missing += 1
            continue
        gff_set = frozenset(grec["alleles"])
        der_set = frozenset((d["a_forward"], d["b_forward"]))
        comp_set = frozenset(complement(b) for b in der_set)
        row = dict(d)
        row["position"] = str(grec["pos"])
        if gff_set == der_set:
            pass
        elif gff_set == comp_set:
            row["a_forward"], row["b_forward"] = complement(d["a_forward"]), complement(d["b_forward"])
            row["notes"] = (row["notes"] + "; " if row["notes"] else "") + \
                "6K gnm2 alleles complemented vs 50K derivation; A/B complemented"
        else:
            k6_allele_mismatch += 1
            continue
        k6.append(row)

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    for name, rows in (("soysnp50k.csv", derived), ("soysnp6k.csv", k6)):
        with open(DATA_DIR / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=OUT_FIELDS)
            w.writeheader()
            w.writerows(rows)

    print("\n---- derivation summary ----")
    print(f"Table S1 records:          {len(s1)}")
    print(f"  unambiguous A/B:         {stats['unambiguous']}")
    print(f"  sequence-walk A/B:       {stats['walk']}")
    print(f"  walk failures:           {stats['walk_fail']}")
    print(f"  non-ACGT flank:          {stats['non_acgt']}")
    if walk_ks:
        import statistics
        print(f"  walk distance k: max={max(walk_ks)} median={statistics.median(walk_ks)}")
    print(f"gnm1 join: match={stats['gff1_match']} "
          f"allele-mismatch(excluded)={stats['gff1_allele_mismatch']} absent(excluded)={stats['gff1_absent']}")
    print(f"gnm2 join: match={stats['gnm2_match']} complemented={stats['gnm2_complement']} "
          f"mismatch(excluded)={stats['gnm2_mismatch']} absent={stats['gnm2_absent']}")
    print(f"50K lookup rows written:   {len(derived)}")
    print(f"6K: inherited={len(k6)} missing-from-50K={k6_missing} allele-mismatch(excluded)={k6_allele_mismatch}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
