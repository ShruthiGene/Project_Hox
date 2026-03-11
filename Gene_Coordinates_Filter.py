#!/usr/bin/env python3
"""
filter_gene_coords.py

Filter a BED file of gene coordinates to only include genes
present in one or more gene lists (supplied as .bed or .xlsx).

Usage:
    python filter_gene_coords.py \\
        --gene-list  genes.bed  [or genes.xlsx] \\
        --coords     all_gene_coords.bed \\
        --output     filtered_coords.bed \\
        [--columns   0 1 2]          # xlsx: 0-based column indices to use (default: all)
        [--sheet     Sheet1]         # xlsx: sheet name or index (default: first sheet)
        [--no-header]                # xlsx: file has no header row
        [--coord-name-col 3]         # coords BED: 0-based col index of gene name (default: 3)
        [--case-insensitive]         # match gene names case-insensitively

Arguments:
    --gene-list   Path to gene list file (.bed or .xlsx).
                  BED  – every non-blank, non-comment token in the file is treated as a gene name.
                  XLSX – each column becomes a separate gene list; all are merged by default.
    --coords      Path to the reference gene-coordinates BED file
                  (format: chr  start  end  gene_name  [other fields …]).
    --output      Path for the filtered output BED file (default: filtered_coords.bed).
    --columns     (xlsx only) Space-separated 0-based column indices to include.
    --sheet       (xlsx only) Sheet name or 0-based index (default: 0).
    --no-header   (xlsx only) Treat the first row as data, not a header.
    --coord-name-col  0-based column index of the gene name in the coords BED (default: 3).
    --case-insensitive  Ignore case when matching gene names.
"""

import argparse
import os
import sys


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def load_genes_from_bed(path: str, case_insensitive: bool) -> set:
    """Read every whitespace-separated token from a BED-style gene-list file."""
    genes = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            for token in line.split():
                genes.add(token.lower() if case_insensitive else token)
    return genes


def load_genes_from_xlsx(path: str,
                         columns: list | None,
                         sheet,
                         has_header: bool,
                         case_insensitive: bool) -> set:
    """Read gene names from an Excel file (one or more columns)."""
    try:
        import openpyxl
    except ImportError:
        sys.exit(
            "ERROR: openpyxl is required to read .xlsx files.\n"
            "Install it with:  pip install openpyxl"
        )

    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)

    # resolve sheet
    if isinstance(sheet, int):
        ws = wb.worksheets[sheet]
    else:
        ws = wb[sheet]

    genes = set()
    for row_idx, row in enumerate(ws.iter_rows(values_only=True)):
        if row_idx == 0 and has_header:
            continue  # skip header
        for col_idx, cell in enumerate(row):
            if columns is not None and col_idx not in columns:
                continue
            if cell is None:
                continue
            token = str(cell).strip()
            if token:
                genes.add(token.lower() if case_insensitive else token)

    wb.close()
    return genes


def filter_coords(coords_path: str,
                  gene_set: set,
                  name_col: int,
                  case_insensitive: bool,
                  output_path: str) -> tuple[int, int]:
    """
    Stream through the coords BED and write matching lines.
    Returns (total_lines, matched_lines).
    """
    total = matched = 0
    with open(coords_path) as fh_in, open(output_path, "w") as fh_out:
        for line in fh_in:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                fh_out.write(line)          # preserve comments / blank lines
                continue
            total += 1
            fields = stripped.split("\t")
            if name_col >= len(fields):
                # fall back to whitespace split
                fields = stripped.split()
            if name_col >= len(fields):
                print(f"WARNING: line has no column {name_col}, skipping: {stripped[:80]}")
                continue
            gene_name = fields[name_col]
            lookup = gene_name.lower() if case_insensitive else gene_name
            if lookup in gene_set:
                fh_out.write(line)
                matched += 1
    return total, matched


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Filter a gene-coordinates BED file to genes in a gene list.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--gene-list", required=True,
                   help="Gene list file (.bed or .xlsx)")
    p.add_argument("--coords", required=True,
                   help="Reference coordinates BED file (chr start end gene_name …)")
    p.add_argument("--output", default="filtered_coords.bed",
                   help="Output BED file (default: filtered_coords.bed)")
    p.add_argument("--columns", nargs="+", type=int, default=None,
                   help="(xlsx) 0-based column indices to read (default: all columns)")
    p.add_argument("--sheet", default=0,
                   help="(xlsx) Sheet name or 0-based index (default: 0)")
    p.add_argument("--no-header", action="store_true",
                   help="(xlsx) First row is data, not a header")
    p.add_argument("--coord-name-col", type=int, default=3,
                   help="0-based column index of gene name in coords BED (default: 3)")
    p.add_argument("--case-insensitive", action="store_true",
                   help="Match gene names case-insensitively")
    return p.parse_args()


def main():
    args = parse_args()

    # --- validate inputs ---
    for label, path in [("--gene-list", args.gene_list), ("--coords", args.coords)]:
        if not os.path.isfile(path):
            sys.exit(f"ERROR: {label} file not found: {path}")

    ext = os.path.splitext(args.gene_list)[1].lower()

    # resolve sheet (may be int or str)
    sheet = args.sheet
    try:
        sheet = int(sheet)
    except (ValueError, TypeError):
        pass  # keep as string

    # --- load gene list ---
    print(f"Loading gene list from: {args.gene_list}")
    if ext == ".xlsx":
        gene_set = load_genes_from_xlsx(
            args.gene_list,
            columns=args.columns,
            sheet=sheet,
            has_header=not args.no_header,
            case_insensitive=args.case_insensitive,
        )
    elif ext in (".bed", ".txt", ".tsv", ".csv", ""):
        gene_set = load_genes_from_bed(args.gene_list, args.case_insensitive)
    else:
        print(f"WARNING: unrecognised extension '{ext}', treating as BED/text.")
        gene_set = load_genes_from_bed(args.gene_list, args.case_insensitive)

    print(f"  → {len(gene_set):,} unique gene names loaded")

    if not gene_set:
        sys.exit("ERROR: No gene names were loaded from the gene-list file. Aborting.")

    # --- filter coords BED ---
    print(f"Filtering coordinates BED: {args.coords}")
    total, matched = filter_coords(
        coords_path=args.coords,
        gene_set=gene_set,
        name_col=args.coord_name_col,
        case_insensitive=args.case_insensitive,
        output_path=args.output,
    )

    print(f"  → {matched:,} / {total:,} coordinate entries matched")
    print(f"Output written to: {args.output}")

    if matched == 0:
        print(
            "\nHINT: No genes matched. Check that:\n"
            f"  • Column {args.coord_name_col} in '{args.coords}' really contains gene names.\n"
            "  • Gene name formatting matches between the two files.\n"
            "  • Try --case-insensitive if capitalisation differs."
        )


if __name__ == "__main__":
    main()
