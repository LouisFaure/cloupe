import argparse
from pprint import pprint
from cloupe import Cloupe


def main():
    parser = argparse.ArgumentParser(prog="Cloupe parser", description="Read Cloupe file")
    parser.add_argument("-c", "--cloupe", required=True, type=str, help="Input Cloupe file")

    parser.add_argument(
        "--csv",
        action="store_true",
        help="Export barcodes, features, observations and projections to CSV",
    )

    parser.add_argument(
        "--anndata",
        action="store_true",
        help="Export cloupe to anndata",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        default=".",
        required=False,
        type=str,
        help="Control the output folder to write the extracted data as CSV",
    )

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        SystemExit()

    c = Cloupe(args.cloupe)
    pprint(c)

    if args.csv:
        c.to_csv(outdir=args.outdir)
    elif args.anndata:
        c.to_anndata()
