#!/usr/bin/env python3
#This script generate a table of alternative transcript identifiers per locus, one line per locus. 
#It is intended to be used to prepare inputs for OMARK
import argparse
import gzip
from collections import defaultdict


def open_file(path):
    """Open normal or gzipped files transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_attributes(attr):
    """Parse GFF3 attribute column."""
    d = {}
    for item in attr.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def main():
    parser = argparse.ArgumentParser(description="Extract mRNA IDs grouped by Parent from a GFF3 file")
    parser.add_argument("input_gff", help="Input GFF3 file (can be .gz)")
    parser.add_argument("output_file", help="Output file")
    args = parser.parse_args()

    parent_dict = defaultdict(set)

    with open_file(args.input_gff) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")
            if len(cols) < 9:
                continue

            if cols[2] != "mRNA" and cols[2] != 'transcript':
                continue

            attrs = parse_attributes(cols[8])
            ID = attrs.get("ID")
            Parent = attrs.get("Parent")

            if ID and Parent:
                parent_dict[Parent].add(ID)

    with open(args.output_file, "w") as out:
        for parent in parent_dict:
            ids = sorted(parent_dict[parent])
            out.write(";".join(ids) + "\n")


if __name__ == "__main__":
    main()
