#!/usr/bin/env python3
import sys
import csv
from pathlib import Path
import gzip
import argparse
import pandas as pd
from utils import read_bed_file

def load_mapping(mapping_file):
    mapping_file = Path(mapping_file)
    # Read the mapping file as a DataFrame
    df = pd.read_csv(mapping_file, sep="\t")
    # Check if the required columns exist
    required_cols = ['UCSC.style.name', 'Sequence.Name']
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Error: Mapping file is missing required column(s): {', '.join(missing)}")
    # Build a mapping from NCBI-style (Sequence.Name) to UCSC-style (UCSC.style.name)
    mapping = dict(zip(df['Sequence.Name'], df['UCSC.style.name']))
    return mapping


def main():
    parser = argparse.ArgumentParser(description="Convert a BED file from NCBI to UCSC style if needed")
    parser.add_argument("--bed", "-b", required=True, help="Input BED file (or gzipped BED file)")
    parser.add_argument("--mapping",  "-m", required=True, help="Mapping file with columns 'UCSC.style.name' and 'Sequence.Name'")
    parser.add_argument("--output", "-o", required=True, help="Output converted BED file")
    
    args = parser.parse_args()
    bed_path = Path(args.bed)
    mapping_path = Path(args.mapping)
    out_path = Path(args.output)

    bed_df = read_bed_file(bed_path)

    # Define recognized sets:
    ucsc_set = { "chr" + str(i) for i in range(1,23) }
    ucsc_set.update({"chrX", "chrY"})
    ncbi_set = { str(i) for i in range(1,23) }
    ncbi_set.update({"X", "Y"})

    # Scan all non-empty lines to see if any chromosome name is in UCSC or NCBI style.
    found_ucsc = False
    found_ncbi = False

    unique_chroms = set(bed_df["chrom"].dropna().unique())
    if any(ch in ucsc_set for ch in unique_chroms):
        found_ucsc = True
    elif any(ch in ncbi_set for ch in unique_chroms):
        found_ncbi = True

    if found_ucsc and not found_ncbi:
        # If only UCSC-style names are found, no conversion is needed.
        needs_conversion = False
    elif found_ncbi and not found_ucsc:
        # If only NCBI-style names are found, conversion is needed.
        needs_conversion = True
    else:
        sys.exit("Error: No recognizable chromosome names found in BED file (neither UCSC nor NCBI style).")

    needs_conversion

    # If conversion is needed, load the mapping and convert the 'chrom' column.
    if needs_conversion:
        mapping = load_mapping(mapping_path)
        def convert_chrom(ch):
            if ch in mapping:
                return mapping[ch]
            else:
                sys.stderr.write(f"Warning: Chromosome '{ch}' not found in mapping file. Leaving it unchanged.\n")
                return ch
        bed_df["chrom"] = bed_df["chrom"].apply(convert_chrom)

    # Write the (possibly converted) BED file.
    bed_df.to_csv(out_path, sep="\t", index=False, header=False)
    return 0

if __name__ == "__main__":
    sys.exit(main())