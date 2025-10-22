#!/usr/bin/env python3
"""
long_to_matrix_strip_prefix.py

Usage:
    python long_to_matrix_strip_prefix.py unfiltered_mash.csv out_matrix.csv

This script expects a header with these columns:
    assembly_file,reference,mash_distance

It will remove the prefix 'bb_cluster_references/' from the reference values before pivoting.
"""

import sys
import pandas as pd

PREFIX = "bb_cluster_references/"

def clean_reference(x, prefix=PREFIX):
    if pd.isna(x):
        return x
    s = str(x)
    if s.startswith(prefix):
        return s[len(prefix):]
    return s

def long_to_matrix_strip_prefix(infile, outfile, fill_value=""):
    df = pd.read_csv(infile, dtype=str)
    required = {"assembly_file", "reference", "mash_distance"}
    if not required.issubset(df.columns):
        raise ValueError(f"Input file must contain columns: {', '.join(required)}")

    # Clean reference values by removing the prefix
    df["reference"] = df["reference"].apply(clean_reference)

    # Convert mash_distance to numeric if possible
    df["mash_distance"] = pd.to_numeric(df["mash_distance"], errors="coerce")

    # Pivot to matrix
    mat = df.pivot_table(index="assembly_file",
                         columns="reference",
                         values="mash_distance",
                         aggfunc="first")

    # Sort rows/columns for stable output
    mat = mat.sort_index(axis=0).sort_index(axis=1)

    # Write to CSV; keep empty cells blank by default
    mat.to_csv(outfile, float_format="%.6f", na_rep="")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python long_to_matrix_strip_prefix.py input_long.csv output_matrix.csv")
        sys.exit(1)
    infile = sys.argv[1]
    outfile = sys.argv[2]
    long_to_matrix_strip_prefix(infile, outfile)
