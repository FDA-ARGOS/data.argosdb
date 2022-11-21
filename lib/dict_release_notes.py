#!/usr/bin/env python3
"""Data Dictionary Release Notes

Generates release notes for ARGOS Data Dictionary.
"""

import csv
from argparse import ArgumentParser, SUPPRESS
import sys
import os
__version__ = "0.1.0"
__status__ = "Production"

def usr_args():
    """User Arguments

    User supplied arguments from command line for function

    Returns
    -------
        ArgumentParser objects to be digested by subsequent functions.
    """

    parser = ArgumentParser(
        add_help=False,
        prog='data_sheet_validator',
        description="Release Note Generator. "
            "Used to generate release notes from two differnt versions.")

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-o', '--old',
        required=True,
        help="Directory for old version of dictionary.")

    required.add_argument('-n', '--new',
        required=True,
        help="Directory for new version of dictionary")

    optional.add_argument('-v', '--version',
        action='version',
        version='%(prog)s ' + __version__)
    optional.add_argument('-h', '--help',
        action='help',
        default=SUPPRESS,
        help='show this help message and exit')

    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    return parser.parse_args()


def load_tsv(options):
    """Load TSVs

    Parameters
    ----------
    options.old: str
        An inpit directory string.
    options.new: str
        An inpit directory string.

     Returns
    -------
        Printout of comparison stats.
    """
    counts ={}
    readme = []
    property_definition = {}
    core_property_list = []
    non_core_property_list = []
    new_readme = []
    new_property_definition = {}
    new_core_property_list = []
    new_non_core_property_list = []

    for sheet in os.listdir(options.old):
        if 'README' in sheet:
            file_path = os.path.join(options.old, sheet)
            with open(file_path, 'r', encoding='utf-8') as file:
                for row in file:
                    if row == '\t\t\n':
                        continue
                    readme.append(row)

        if 'property_definition' in sheet:
            file_path = os.path.join(options.old, sheet)
            with open(file_path, 'r', encoding='utf-8') as file:
                data = csv.reader(file, delimiter='\t')
                next(data)
                for row in data:
                    property_definition[row[0]] = row[1:]

        if 'core_property_list' in sheet:
            file_path = os.path.join(options.old, sheet)
            with open(file_path, 'r', encoding='utf-8') as file:
                data = csv.reader(file, delimiter='\t')
                next(data)
                counts['core_property_list'] ={}
                for row in data:
                    core_property_list.append(row)
                    if row[1] not in counts['core_property_list']:
                        counts['core_property_list'][row[1]] = 1
                    else:
                        counts['core_property_list'][row[1]] += 1

        if 'annotation_property_list' in sheet or\
            'non_core_property_list' in sheet:
            file_path = os.path.join(options.old, sheet)
            with open(file_path, 'r', encoding='utf-8') as file:
                data = csv.reader(file, delimiter='\t')
                next(data)
                counts['non_core_property_list'] ={}
                for row in data:
                    non_core_property_list.append(row)
                    if row[1] not in counts['non_core_property_list']:
                        counts['non_core_property_list'][row[1]] = 1
                    else:
                        counts['non_core_property_list'][row[1]] += 1

    for sheet in os.listdir(options.new):
        if 'README' in sheet:

            file_path = os.path.join(options.new, sheet)
            with open(file_path, 'r', encoding='utf-8') as file:
                for row in file:
                    if row == '\t\t\n':
                        continue
                    new_readme.append(row)

        if 'property_definition.tsv' in sheet:
            file_path = os.path.join(options.new, sheet)
            with open(file_path, 'r', encoding='utf-8') as file:
                data = csv.reader(file, delimiter='\t')
                next(data)
                for row in data:
                    new_property_definition[row[0]] = row[1:]

        if 'core_property_list' in sheet:
            file_path = os.path.join(options.new, sheet)
            with open(file_path, 'r', encoding='utf-8') as file:
                data = csv.reader(file, delimiter='\t')
                next(data)
                counts['new_core_property_list'] ={}
                for row in data:
                    new_core_property_list.append(row)
                    if row[1] not in counts['new_core_property_list']:
                        counts['new_core_property_list'][row[1]] = 1
                    else:
                        counts['new_core_property_list'][row[1]] += 1

        if 'annotation_property_list' in sheet or\
            'non-core_property_list' in sheet:
            file_path = os.path.join(options.new, sheet)
            with open(file_path, 'r', encoding='utf-8') as file:
                data = csv.reader(file, delimiter='\t')
                next(data)
                counts['new_non_core_property_list'] ={}
                for row in data:
                    new_non_core_property_list.append(row)
                    if row[1] not in counts['new_non_core_property_list']:
                        counts['new_non_core_property_list'][row[1]] = 1
                    else:
                        counts['new_non_core_property_list'][row[1]] += 1
    print(counts)
    removed = []
    added = []
    # print(new_non_core_property_list == non_core_property_list)
    # print(readme == new_readme)
    # print(property_definition == new_property_definition)
    # print(core_property_list == new_core_property_list)
    # print(non_core_property_list == new_non_core_property_list)
    # print('\n')
    for key in property_definition:
        if key not in new_property_definition.keys():
            removed.append(key)
            # print(f'{key} removed from {options.new}')

    for key in new_property_definition:
        if key not in property_definition.keys():
            added.append(key)
            # print(f'{key} added to {options.new}')

    print('\nCore Property Sheet Stats:\n Sheet \t Properties')
    for key in counts['new_core_property_list']:
        print(key, '\t', counts['new_core_property_list'][key] )

    print('\nNon Core Property Sheet Stats:\n Sheet \t Properties')
    for key in counts['new_non_core_property_list']:
        print(key, '\t', counts['new_non_core_property_list'][key] )
    print('\n')
    print(f'{len(removed)} items removed from the dictionary: {removed}')
    print('\n')
    print(f'{len(added)} items added to the dictionary: {added}')

def main():
    """
    Main function
    """
    options = usr_args()
    load_tsv(options)

#______________________________________________________________________________#
if __name__ == "__main__":
    main()
