"""Validate the bundled BARCSoySNP6K lookup table against an official manifest.

For every SNP present in both the Illumina manifest and ``soysnp6k.csv``
(joined by SNP name), this checks that the manifest's ``SNP`` column
``[X/Y]`` -- allele A / allele B on the assay's ``IlmnStrand`` -- is
consistent with the table's ``(a_forward, b_forward)``:

    (X, Y) == (a_forward, b_forward)  or  (X, Y) == complement(a_forward, b_forward)

An A/B label swap (``(X, Y) == (b_forward, a_forward)``) counts as a
disagreement.  The check is strand-aware but label-strict: a SNP passes only
if the manifest's A allele is the table's A allele (possibly complemented),
never the table's B allele.

Usage:
    python validate_against_manifest.py MANIFEST.csv [--table soysnp6k.csv]

Exits 0 when every joined SNP agrees, 1 otherwise.
"""

import argparse
import csv
import re
import sys

COMP = str.maketrans("ACGT", "TGCA")


def read_manifest(path):
    """Return {snp_name: (A_allele, B_allele, ilmn_strand, customer_strand)}."""
    snps = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        in_assay, header = False, None
        for row in reader:
            if not row:
                continue
            if row[0] == "[Assay]":
                in_assay = True
                continue
            if not in_assay:
                continue
            if row[0] == "IlmnID":
                header = row
                continue
            d = dict(zip(header, row))
            m = re.match(r"\[([ACGT])/([ACGT])\]", d["SNP"])
            if not m:
                print(f"warning: unparseable SNP column {d['SNP']!r} "
                      f"for {d['IlmnID']}", file=sys.stderr)
                continue
            snps["BARC_1.01_" + d["Name"]] = (
                m.group(1), m.group(2), d["IlmnStrand"], d["CustomerStrand"])
    return snps


def read_table(path):
    """Return {snp_id: (a_forward, b_forward)}."""
    table = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            table[row["snp"]] = (row["a_forward"], row["b_forward"])
    return table


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("manifest", help="Illumina manifest CSV ([Assay] section)")
    ap.add_argument("--table", default="soysnp6k.csv",
                    help="bundled lookup table CSV")
    args = ap.parse_args()

    manifest = read_manifest(args.manifest)
    table = read_table(args.table)
    joined = [k for k in manifest if k in table]
    manifest_only = [k for k in manifest if k not in table]
    table_only = [k for k in table if k not in manifest]

    agree, disagreements = 0, []
    for key in joined:
        x, y, ilmn, _cust = manifest[key]
        af, bf = table[key]
        if (x, y) == (af, bf) or (x, y) == (af.translate(COMP),
                                            bf.translate(COMP)):
            agree += 1
        else:
            disagreements.append((key, x, y, af, bf, ilmn))

    print(f"manifest SNPs: {len(manifest)}")
    print(f"table SNPs:    {len(table)}")
    print(f"joined:        {len(joined)}")
    print(f"agree:         {agree}")
    print(f"disagree:      {len(disagreements)}")
    print(f"manifest-only: {len(manifest_only)}")
    print(f"table-only:    {len(table_only)}")
    for key, x, y, af, bf, ilmn in disagreements[:20]:
        print(f"  DISAGREE {key}: manifest [{x}/{y}] "
              f"vs table ({af},{bf}) IlmnStrand={ilmn}")
    if manifest_only:
        print("  manifest-only (first 10):", manifest_only[:10])

    return 0 if not disagreements else 1


if __name__ == "__main__":
    sys.exit(main())
