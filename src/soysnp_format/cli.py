"""Command-line interface for soysnp-format."""

import argparse
import json
import sys

from . import formats
from .convert import convert_genotypes, write_long_table, write_wide_matrix, ConversionError
from .detect import detect_format
from .io import read_genotypes
from .lookup import BUNDLED_CHIPS, bundled_snp_sets, load_lookup


def _add_input_args(p):
    p.add_argument("input", help="genotype file (GenomeStudio export or SNP matrix)")
    p.add_argument("--layout", default="auto",
                   choices=["auto", "final-report", "wide"],
                   help="input layout (default: auto-detect)")
    p.add_argument("--chip", default="auto", choices=["auto"] + list(BUNDLED_CHIPS),
                   help="bundled chip table for SNP matching (default: auto = both)")
    p.add_argument("--lookup",
                   help="custom lookup CSV (snp,a_forward,b_forward) instead of bundled tables")


def _load_lookup(args):
    if args.lookup:
        return load_lookup(path=args.lookup)
    try:
        return load_lookup(args.chip)
    except (FileNotFoundError, ValueError) as exc:
        print(f"warning: {exc}", file=sys.stderr)
        return None


def cmd_detect(args):
    snp_sets = ()
    probe = _load_lookup(args)
    if probe is not None:
        snp_sets = (probe.snps(),)
    data = read_genotypes(args.input, layout=args.layout, snp_sets=snp_sets)
    result = detect_format(data, lookup=probe)
    if args.json:
        print(json.dumps(result.as_dict(), indent=2))
    else:
        print(result.summary())
    return 0


def cmd_convert(args):
    snp_sets = ()
    lookup = _load_lookup(args)
    if lookup is not None:
        snp_sets = (lookup.snps(),)
    data = read_genotypes(args.input, layout=args.layout, snp_sets=snp_sets)
    try:
        out = convert_genotypes(
            data,
            from_format=args.from_.upper(),
            to_format=args.to.upper(),
            lookup=lookup,
            on_unknown_snp=args.on_unknown_snp,
            on_bad_allele=args.on_bad_allele,
        )
    except ConversionError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if args.out_layout == "long":
        enc_name = "AB" if args.to.upper() == "AB" else "Forward"
        write_long_table(out, args.output, delimiter=args.delimiter, encoding_name=enc_name)
    else:
        write_wide_matrix(out, args.output, delimiter=args.delimiter)
    n_snps = len(out.snps)
    print(f"wrote {args.output}: {n_snps} SNPs x {len(out.samples)} samples "
          f"({args.from_.upper()} -> {args.to.upper()})")
    return 0


def cmd_chips(_args):
    print("bundled chip tables:")
    for chip in BUNDLED_CHIPS:
        try:
            lookup = load_lookup(chip)
            print(f"  {chip}: {len(lookup)} SNPs")
        except FileNotFoundError as exc:
            print(f"  {chip}: MISSING ({exc})")
    return 0


def build_parser():
    p = argparse.ArgumentParser(
        prog="soysnp-format",
        description="Detect and translate allele encodings (Illumina A/B vs forward "
                    "strand) in soybean SNP array exports.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    d = sub.add_parser("detect", help="report the likely allele encoding of a file")
    _add_input_args(d)
    d.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    d.set_defaults(func=cmd_detect)

    c = sub.add_parser("convert", help="translate a file between allele encodings")
    _add_input_args(c)
    c.add_argument("--from", dest="from_", default="auto",
                   help="source encoding: AB, FORWARD, or auto (default: auto-detect)")
    c.add_argument("--to", required=True,
                   help="target encoding: AB or FORWARD")
    c.add_argument("-o", "--output", required=True, help="output file path")
    c.add_argument("--out-layout", default="wide", choices=["wide", "long"],
                   help="output layout (default: wide SNP-rows matrix)")
    c.add_argument("--delimiter", default="\t", help="output delimiter (default: tab)")
    c.add_argument("--on-unknown-snp", default="error", choices=["error", "skip"],
                   help="SNPs absent from the lookup table (default: error)")
    c.add_argument("--on-bad-allele", default="error", choices=["error", "missing"],
                   help="genotypes outside the source alphabet (default: error)")
    c.set_defaults(func=cmd_convert)

    k = sub.add_parser("chips", help="list bundled chip lookup tables")
    k.set_defaults(func=cmd_chips)

    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
