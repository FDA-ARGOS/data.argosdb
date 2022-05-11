#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Header Parser

usage: argosdb [options]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -c CSV, --csv CSV     CSV file to create.
  -s SCHEMA, --schema SCHEMA
                        Root json schema to parse
"""

import csv
import json
import argparse
import jsonschema
from jsonschema import validate
import sys

__version__ = "1.1.0"
__status__ = "Production"

def usr_args():
    """User Arguments

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

    parser.add_argument('-t', '--txt',
                                required=True,
                                help="Text file to read.")

    parser.add_argument('-d', '--dict',
                                # type = argparse.FileType('r'),
                                help="Property list to modify")

    # Print usage message if no args are supplied.
    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    options = parser.parse_args()
    return options


def read_headers(options):
    """Read Headers

    Reads a header file output by `head -1 * > ~/headers.txt` at
    `argosdb-vm-dev:/data/shared/argosdb/generated/datasets/reviewed/`
    """

    files = {}
    file_names = []
    property_file_names = {}
    with open(options.txt, 'r', encoding='utf-8') as file:
        # lines = file.readlines()
        for line in file:
            if line.startswith('==> '):
                line = line.replace('==> ', '')
                line = line.replace(' <==', '')
                line = line.replace('\n', '')
                file_names.append(line)
                headerline = next(file)
                if headerline.startswith('>'):
                    continue
                if headerline.startswith('"'):
                    headerline = headerline.replace('"','')
                    headerline = headerline.replace('\n', '')
                    files[line] = headerline.split(',')
                    for header in headerline.split(','):
                        if header in property_file_names.keys():
                            property_file_names[header].append(line)
                        else:
                            property_file_names[header] = [line]
                else:
                    headerline = headerline.replace('\n', '')
                    files[line] = headerline.split('\t')
                    for header in headerline.split('\t'):
                        if header in property_file_names.keys():
                            property_file_names[header].append(line)
                        else:
                            property_file_names[header] = [line]
                    # print(headerline)

        # print(json.dumps(property_file_names))
        # print(len(property_file_names.keys()))
        # print(len(file_names))
        return file_names, property_file_names, files

def read_property_list(options, file_names, property_file_names, files):
    """
        This function writes the schema from the given csv
    """
    # print(options.dict)
    with open(options.dict, 'r', encoding='utf-8') as file:
        data = csv.reader(file, delimiter="\t")
        header = next(data)
        for row in data:
            file_list = ''
            try:
                for prop in property_file_names[row[0]]:
                    file_list += str(prop + '|')
                file_list = file_list.removesuffix('|')
                print(row[0], '\t', file_list)
            except KeyError:
                print(f'ERROR! {row[0]}')

def main():
    """
    Main function
    """

    options = usr_args()
    file_names, property_file_names, files = read_headers(options)
    read_property_list(options, file_names, property_file_names, files)

if __name__ == "__main__":
    main()
