#!/usr/bin/env python3
"""ARGOS Schema Generator

This script will perform various functions related to processing the ARGOS
data dictionary and associated schemas. General help below.

usage: write_schema.py -i INPUT [-o OUTPUT] [-v] [-h]

Schema Creation. Used to create a schema from a data sheet.

required arguments:
  -i INPUT, --input INPUT
                        Directory for schema files. This directory should contain 'core_property_list.tsv','non_core_property_list.tsv', and 'property_definition.tsv'.These files can be a
                        'tsv' or 'csv'. Other file types are not permitted.

optional arguments:
  -o OUTPUT, --output OUTPUT
                        Output file to create. Default is a JSON file.
  -v, --version         show program's version number and exit
  -h, --help            show this help message and exit
"""

import csv
import json
from argparse import ArgumentParser, SUPPRESS
import sys
import os

__version__ = "0.8"
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
        prog='write_schema.py',
        description="Schema Creation. "
            "Used to create a schema from a data sheet. ")

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-i', '--input',
        required=True,
        help="Directory for schema files. "
             "This directory should contain 'core_property_list.tsv',"
             "'non_core_property_list.tsv', and 'property_definition.tsv'."
             "These files can be a 'tsv' or 'csv'. "
            "Other file types are not permitted." )

    optional.add_argument('-o', '--output',
        help="Output file to create. Default is a JSON file. ")
    # optional.add_argument('-m', '--multi',
    #     action='store_true',
    #     help='Flag to indicate if multiple items are being processed')
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

def list_2_schema(options):
    """Create Schema JSON

    Takes a TSV and writes JSON formatted schema

    Parameters
    ----------
    options.input: str
        An inpit file to create the schema/schemas. This should be a TSV.

    definition: str
        a list of property definitions.

    options.output: str, optional
        An output file. If this is supplied then the function output will be
        written to this directory.

    Returns
    -------
        Either writes JSON schema object/objects to a file (if options.output
        is supplied) or prints to the terminal.
    """

    core_sheet = options.input + 'core_property_list.tsv'
    non_core_sheet = options.input + 'non-core_property_list.tsv'
    definitions = options.input + 'property_definition.tsv'
    output_dir = 'schema/v'+__version__+'/'
    write_output = True
    raw_url = 'https://raw.githubusercontent.com/FDA-ARGOS/data.argosdb/v'\
            +__version__+'/schema/'+__version__+'/'
    if not os.path.exists(core_sheet):
        print(f'ERROR!! \'{core_sheet}\' does not exist')
        exit()
    if not os.path.exists(non_core_sheet):
        print(f'ERROR!! \'{non_core_sheet}\' does not exist')
        exit()
    if not os.path.exists(definitions):
        print(f'ERROR!! \'{definitions}\' does not exist')
        exit()

    if not options.output:
        print("No output directory supplied [-o]. Printing to terminal")
        write_output = False
    prop_defs = {}

    with open(definitions, 'r', encoding='utf8') as definitions:
        def_data = csv.reader(definitions, delimiter="\t")
        for row in def_data:
            prop_defs[row[0].rstrip()] = row[3]
    print(write_output, raw_url)

    with open(core_sheet, 'r', encoding='utf8') as file:
        argos_schemas = {}
        data = csv.reader(file, delimiter="\t")
        next(data)
        count = 0
        for row in data:
            count+= 1
            try:
                prop_defs[row[0].rstrip()]
            except KeyError:
                print(f"Error! {row[0]} at row {count} in {core_sheet} is not defined in\
                    {definitions.name}. Exiting.")
                return
            if row[1] not in argos_schemas:
                argos_schemas[row[1]] = {
                    'definitions': {},
                    '$schema': 'http://json-schema.org/draft-07/schema#',
                    '$id': raw_url+row[1],
                    'title': row[1],
                    'type': 'object',
                    'required':[],
                    'properties':{}
                }
            argos_schemas[row[1]]['properties'][row[0]] = {
                '$id':row[3],
                'title': row[4],
                'description': prop_defs[row[0].rstrip()],
                'type': row[5],
                'default': row[6],
                'examples': [row[7]],
                'pattern': row[8]
            }
            if row[2] == 'required':
                (argos_schemas[row[1]]['required']).append(row[0])
        jsonf = json.dumps(argos_schemas, indent=4)
        for item in argos_schemas:
            file_name = output_dir+'core/'+item.split('.')[0]+'.json'
            with open(file_name, 'w', encoding='utf-8') as file:
                file.write(json.dumps(argos_schemas[item], indent=4))

    with open(non_core_sheet, 'r', encoding='utf8') as file:
        argos_schemas = {}
        data = csv.reader(file, delimiter="\t")
        next(data)
        count = 1
        for row in data:
            count+= 1
            try:
                prop_defs[row[0].rstrip()]
            except KeyError:
                print(f"Error! {row[0]} at row {count} in {non_core_sheet} is not defined in \
                    {definitions.name}. Exiting.")
                return
            if row[1] not in argos_schemas:
                argos_schemas[row[1]] = {
                    'definitions': {},
                    '$schema': 'http://json-schema.org/draft-07/schema#',
                    '$id': raw_url+row[1],
                    'title': row[1],
                    'type': 'object',
                    'required':[],
                    'properties':{}
                }
            argos_schemas[row[1]]['properties'][row[0]] = {
                '$id':row[3],
                'title': row[4],
                'description': prop_defs[row[0].rstrip()],
                'type': row[5],
                'default': row[6],
                'examples': [row[7]],
                'pattern': row[8]
            }
            if row[2] == 'required':
                (argos_schemas[row[1]]['required']).append(row[0])
        jsonf = json.dumps(argos_schemas, indent=4)
        for item in argos_schemas:
            file_name = output_dir+'non-core/'+item.split('.')[0]+'.json'
            with open(file_name, 'w', encoding='utf-8') as file:
                file.write(json.dumps(argos_schemas[item], indent=4))

def main():
    """
    Main function
    """

    options = usr_args()
    list_2_schema(options)

if __name__ == "__main__":
    main()
