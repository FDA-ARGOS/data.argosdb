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
import sys
import argparse

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
        required=True,
        # type = argparse.FileType('r'),
        help="Property list to modify")

    # Print usage message if no args are supplied.
    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    options = parser.parse_args()
    return options


def read_headers(filepath):
    """Read Headers

    Reads a header file output by `head -1 * > ~/headers.txt` at
    `argosdb-vm-dev:/data/shared/argosdb/generated/datasets/reviewed/`

    Parameters
    ----------
    filepath : str
        A file path for the headers file to be processed

    Returns
    -------
    file_names : list
        a list of file names
    property_file_names : dict
        a dict
    file_headers : dict
        a dict
    """

    file_headers = {}
    file_names = []
    property_file_names = {}
    with open(filepath, 'r', encoding='utf-8') as file:
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
                    file_headers[line] = headerline.split(',')
                    for property_name in file_headers[line]:
                        if property_name in property_file_names:
                            property_file_names[property_name].append(line)
                        else:
                            property_file_names[property_name] = [line]
                else:
                    headerline = headerline.replace('\n', '')
                    file_headers[line] = headerline.split('\t')
                    for property_name in headerline.split('\t'):
                        if property_name in property_file_names:
                            property_file_names[property_name].append(line)
                        else:
                            property_file_names[property_name] = [line]

        return file_names, property_file_names, file_headers

def read_property_list(filepath, property_file_names):
    """Read Property List

    This function writes the schema from the given csv

    Parameters
    ----------
    filepath : str
        A file path for the headers file to be processed
    property_file_names : dict
        a dictionary with properties as the key and a list of the file the
        property is present in as the value

    Returns
    -------

    """

    with open(filepath, 'r', encoding='utf-8') as file:
        data = csv.reader(file, delimiter="\t")
        header = next(data)
        properties = []
        for row in data:
            file_list = ''
            properties.append(row[0])
            try:
                for prop in property_file_names[row[0]]:
                    file_list += str(prop + '|')
                file_list = file_list.removesuffix('|')
                print(row[0], '\t', file_list)
            except KeyError:
                print(row[0], '\t', row[1], '\t', 'ERROR! not found in data sheets')
        print('\n')
        for key in property_file_names:
            if key not in properties:
                file_list = ''
                for file in property_file_names[key]:
                    file_list += str(file + '|')
                file_list = file_list.removesuffix('|')
                print(key, '\t', file_list)

def make_property_list(file_names, file_headers):
    """Make Property List

    """

    print('\n')
    for key in file_headers:
        # print(key)
        for column in file_headers[key]:
            print(column, '\t', key)
    # print(json.dumps(file_names))
    # print(json.dumps(file_headers))
    print('All Files: ', len(file_names))
    print('Data Sheets: ', len(file_headers))

def main():
    """
    Main function
    """

    options = usr_args()
    file_names, property_file_names, file_headers = read_headers(options.txt)
    read_property_list(options.dict, property_file_names)
    make_property_list(file_names, file_headers)

if __name__ == "__main__":
    main()
