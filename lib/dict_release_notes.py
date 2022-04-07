#!/usr/bin/env python3
"""Data Dictionary Release Notes

Generates release notes for ARGOS Data Dictionary.
"""

import csv
import argparse
import sys
import os
__version__ = "0.1.0"
__status__ = "Production"

def usr_args():
    """
    functional arguments for process
    https://stackoverflow.com/questions/27529610/call-function-based-on-argparse
    """

    parser = argparse.ArgumentParser()

    # set usages options
    parser = argparse.ArgumentParser(
        prog='argosdb',
        usage='%(prog)s [options]')

    # version
    parser.add_argument(
                        '-v', '--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    parser.add_argument('-o', '--old',
                        required=True,
                        help="Directory for old version of dictionary.")

    parser.add_argument('-n', '--new',
                        required=True,
                        help="Directory for new version of dictionary")

    # Print usage message if no args are supplied.
    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    options = parser.parse_args()
    return options


def load_tsv(options):
    """
        load tsvs
    """
    counts ={}
    readme = []
    property_definitions = {}
    core_property_list = []
    non_core_property_list = []
    new_readme = []
    new_property_definitions = {}
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

        if 'property_definitions' in sheet:
            file_path = os.path.join(options.old, sheet)
            with open(file_path, 'r', encoding='utf-8') as file:
                data = csv.reader(file, delimiter='\t')
                property_definitions_header = next(data)
                for row in data:
                    property_definitions[row[0]] = row[1:]

        if 'core_property_list' in sheet:
            file_path = os.path.join(options.old, sheet)
            with open(file_path, 'r', encoding='utf-8') as file:
                data = csv.reader(file, delimiter='\t')
                core_property_list_header = next(data)
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
                non_core_property_list_header = next(data)
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

        if 'property_definitions.tsv' in sheet:
            file_path = os.path.join(options.new, sheet)
            with open(file_path, 'r', encoding='utf-8') as file:
                data = csv.reader(file, delimiter='\t')
                new_property_definitions_header = next(data)
                for row in data:
                    new_property_definitions[row[0]] = row[1:]

        if 'core_property_list' in sheet:
            file_path = os.path.join(options.new, sheet)
            with open(file_path, 'r', encoding='utf-8') as file:
                data = csv.reader(file, delimiter='\t')
                new_core_property_list_header = next(data)
                counts['new_core_property_list'] ={}
                for row in data:
                    new_core_property_list.append(row)
                    if row[1] not in counts['new_core_property_list']:
                        counts['new_core_property_list'][row[1]] = 1
                    else:
                        counts['new_core_property_list'][row[1]] += 1

        if 'annotation_property_list' in sheet or\
            'non_core_property_list' in sheet:
            file_path = os.path.join(options.new, sheet)
            with open(file_path, 'r', encoding='utf-8') as file:
                data = csv.reader(file, delimiter='\t')
                new_non_core_property_list_header = next(data)
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
    # print(property_definitions == new_property_definitions)
    # print(core_property_list == new_core_property_list)
    # print(non_core_property_list == new_non_core_property_list)
    # print('\n')
    for key in property_definitions:
        if key not in new_property_definitions.keys():
            removed.append(key)
            # print(f'{key} removed from {options.new}')

    for key in new_property_definitions:
        if key not in property_definitions.keys():
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
    # make_schema(options, schema)



#______________________________________________________________________________#
if __name__ == "__main__":
    main()
