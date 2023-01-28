#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Add Lineage

This script will take in a TSV, a cell number for taxonomy ID, and a cell
number for the lineage, and then write the taxonomic lineages in the file.
"""

import os
import sys
import csv
import argparse

__version__ = "0.1"
__status__ = "Draft"

def usr_args() -> argparse.ArgumentParser:
    """Program Arguments

    User supplied arguments from command line for function

    Returns
    -------
        ArgumentParser objects to be digested by subsequent functions.
    """

    parser = argparse.ArgumentParser(
        add_help=False,
        prog='add_lineage.py',
        usage='%(prog)s [options]')

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-d', '--data',
        type=str,
        required=True,
        help="Data sheet to read and write.")
    required.add_argument('-m', '--mapping',
        type=str,
        required=True,
        help="Mapping sheet to read.")
    required.add_argument('-t', '--taxonomy',
        type=int,
        required=True,
        help="Column number which has the Taxonomy Id.")
    required.add_argument('-l', '--lineage',
        type=int,
        required=True,
        help="Column number which this program will write the lineage to.")

    optional.add_argument('-o', '--output',
        help="Output file to create. Default is a TSV to currnet directory. ")
    optional.add_argument('-v', '--version',
        action='version',
        version='%(prog)s ' + __version__)
    optional.add_argument('-h', '--help',
        action='help',
        default=argparse.SUPPRESS,
        help='show this help message and exit')

    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    return parser.parse_args()

def load_lineage(mapping: str) -> dict:
    """Load Lineage"""
    lineage_map = dict()
    with open(mapping, 'r', encoding='utf-8') as mapping_file:
        data = csv.reader(mapping_file, delimiter='\t')
        for line in data:
            lineage_map[line[0]] = line[1]
    return lineage_map

def write_lineage (lineage_map: dict, data: str, taxonomy: int, lineage:int):
    """Write Lineage"""
    print(data, taxonomy, lineage)
    tmp_file = "tmp.tsv"
    with open(data, "r", encoding='utf-8') as file, open(
        tmp_file, "w", encoding='utf-8') as out_file:

        reader = csv.reader(file, delimiter='\t')
        writer = csv.writer(out_file, delimiter='\t')
        header = next(reader)
        writer.writerow(header)
        for row in reader:
            row[lineage] = lineage_map[row[taxonomy]]
            writer.writerow(row)
        os.rename(tmp_file, data)

    return 0

def main():
    """
    Main function
    """
    options = usr_args()
    lineage_map = load_lineage(options.mapping)
    write_lineage(lineage_map, options.data, options.taxonomy, options.lineage)

if __name__ == "__main__":
    main()
