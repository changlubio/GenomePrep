from __future__ import annotations

import argparse
import sys


def _process_cmd(args: argparse.Namespace) -> int:
    import os
    from genomeprep import process

    if not os.path.exists(args.infile):
        print(f"Error: file not found: {args.infile}", file=sys.stderr)
        return 1

    result = process(
        infile=args.infile,
        datadir=args.datadir,
        outdir=args.outdir,
        outindex=args.outindex,
    )

    if result.status == "Processed":
        print(f"Success. Output in {args.outdir}/")
        if isinstance(result.info, dict):
            for k, v in result.info.items():
                print(f"  {k}: {v}")
    else:
        print(f"Status: {result.status}")
        print(f"Info:   {result.info}")
        return 1
    return 0


def _download_cmd(args: argparse.Namespace) -> int:
    from genomeprep.data import download

    dest = download(dest=args.dest)
    print(f"Done. Files in: {dest}")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="genomeprep",
        description="Preprocess and quality-control consumer DTC genome files.",
    )
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")

    p_proc = sub.add_parser("process", help="Process a DTC genome file.")
    p_proc.add_argument("infile", help="Input file (raw DTC genome or ZIP)")
    p_proc.add_argument("-d", "--datadir", default=None, help="Path to reference data directory")
    p_proc.add_argument("-o", "--outdir", default="./outputs", help="Output directory")
    p_proc.add_argument("-i", "--outindex", default="outindex", help="Output file index/name")

    p_dl = sub.add_parser("download-data", help="Download reference data files.")
    p_dl.add_argument("--dest", default=None, help="Destination directory (default: ~/.cache/genomeprep)")

    args = parser.parse_args()

    if args.command == "process":
        sys.exit(_process_cmd(args))
    elif args.command == "download-data":
        sys.exit(_download_cmd(args))
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
