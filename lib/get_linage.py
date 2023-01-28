#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Get Lineage

This script will take in a TSV and cell number and and return the
unique taxonomic lineages in the file.
"""

import sys
import csv
import urllib
import argparse
from Bio import Entrez

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
        prog='lineage.py',
        usage='%(prog)s [options]')

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-i', '--input',
        required=True,
        help="Data sheet to read.")
    required.add_argument('-c', '--column',
        type=int,
        help="Column number which has the Taxonomy Id.")

    optional.add_argument('-o', '--output',
        help="Output file to create. Default is output to terminal. ")
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

Entrez.email = "hadley_king@gwu.edu"

def load_taxonomy(input_file: str, column: int) -> list:
    """Load Taxonomy"""
    tax_list = []
    with open(input_file, 'r', encoding='utf-8') as in_file:
        data = csv.reader(in_file, delimiter="\t")
        next(data)
        for line in data:
            tax_id = line[column]
            if tax_id not in tax_list:
                tax_list.append(tax_id)
    return tax_list

def get_lineage(tax_list: list, output:bool=False)-> dict:
    """Get Lineage"""
    lineage = []
    if output is True:
        for tax in tax_list:
            print(tax)
            handle = Entrez.efetch(db='taxonomy', id=tax)
            record2 = Entrez.read(handle)
            lineage.append((tax, record2[0]['Lineage']))
        return lineage

    for tax in tax_list:
        try:
            handle = Entrez.efetch(db='taxonomy', id=tax)
            record2 = Entrez.read(handle)
            print(tax, '\t', record2[0]['Lineage'])
        except urllib.error.HTTPError:
            print(tax, '\t', 'HTTP Error 400: Bad Request')

def write_lineage(lineage: list, outfile: str):
    """Write Lineage
    """
    print(outfile)
    with open(outfile, 'w', encoding='utf-8') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        for taxonomy in lineage:
            writer.writerow([taxonomy[0], taxonomy[1]])

def main():
    """
    Main function
    """
    options = usr_args()
    tax_list = load_taxonomy(options.input, options.column)
    if options.output:
        lineage = get_lineage(tax_list, output=True)
        write_lineage(lineage, options.output)
    else:
        lineage = get_lineage(tax_list)

if __name__ == "__main__":
    main()
