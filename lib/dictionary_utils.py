#!/usr/bin/env python3
"""ARGOS dictionary utilities

This script will perform various functions related to processing the ARGOS
data dictionary and associated schemas. General help below.

positional arguments:
{functions,validate,tsv2json}
    functions           List of all available functions
    validate            Validation options. Used to test a data sheet against a JSON schema. If no schema is supplied will throw an error
    tsv2json            Used to convert a TSV into a JSNO schema. If no mapping file is provided, performs default conversions.

optional arguments:
-h, --help            show this help message and exit
-v, --version         show program's version number and exit
"""

import json, csv
import argparse
import sys

__version__ = "0.1.0"
__status__ = "Production"

def usr_args():
    """User supplied arguments for functions

    """

    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser(
        prog='argosdb_dict_utils',
        usage='%(prog)s [options]')

    # version
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s ' + __version__)

    subparsers = parser.add_subparsers()

    # Create parent subparser. Note `add_help=False` & creation via `argparse.`
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-i', '--input',
        required=True,
        help="Input file to process")
    parent_parser.add_argument('-o', '--output',
        help="Output file to create")
    parent_parser.add_argument('-s', '--schema',
        # type = argparse.FileType('r'),
        help="Root json schema to validate against.")
    parent_parser.add_argument('-m', '--multi',
        action='store_true',
        help='Flag to indicate if multiple items are being processed')

    # Create functions subcommand
    parser_list_functions = subparsers.add_parser('functions',
        help="List of all available functions")
    parser_list_functions.set_defaults(func=list_functions)

    # Create a validate subcommand
    parser_validate  = subparsers.add_parser('validate',
        parents=[parent_parser],
        help="Validation options. "
            "Used to test a data sheet against a JSON schema. "
            "If no schema is supplied will throw an error")
    parser_validate.set_defaults(func=validate_schema)

    parser_tsv2json = subparsers.add_parser('tsv2json',
        parents=[parent_parser],
        help="Used to convert a TSV into a JSNO schema."
            " If no mapping file is provided, performs default conversions.")
    parser_tsv2json.set_defaults(func=tsv2json)

    # Print usage message if no args are supplied.
    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    options = parser.parse_args()
    if parser.parse_args().func is list_functions:
        options.func(parser)
    else:
        options.func(options)

def list_functions(parser):
    """List Available Functions

    Parameters
    ----------
    parser: argparse.ArgumentParser
        Argument parser object.

    Returns
    -------
        Prints help for each function
    """

    print('Function List', type(parser))
    subparsers_actions = [
        # pylint: disable=protected-access
        action for action in parser._actions
        # pylint: disable=W0212
        if isinstance(action, argparse._SubParsersAction)]
    # there will probably only be one subparser_action,
    # but better safe than sorry
    for subparsers_action in subparsers_actions:
        # get all subparsers and print help
        for choice, subparser in subparsers_action.choices.items():
            print("Function: '{}'".format(choice))
            print(subparser.format_help())
    print(parser.format_help())

def validate_schema(options):
    """Checks for Schema Complience[WIP]

    Parameters
    ----------
    options.input: str
        Argument parser object holding attributes to process.
        This should inclued an inpit file and may include an optional output
        file and/or an optional schema

    options.schema: str, optional
    
    options.output: str, optional

    Returns
    -------
        Prints help for each function


    """
    print('validate', type(options.input))

def tsv2json(options):
    """Create Schema JSON

    Tales a TSV and writes JSON formatted schema

    Parameters
    ----------
    options.input: str
        An inpit file to create the schema/schemas. This should be a TSV. 

    options.multi: bool
        Default is False. If true the input file is treated as a flat version
        of a multiple schemas. These will output to a list of JSON.

    options.output: str, optional
        An output file. If this is supplied then the function output will be
        written to this file.

    Returns
    -------
        Either writes JSON schema object/objects to a file (if options.output
        is supplied) or prints to the terminal. 
    """

    definition = 'data_files/property_definition.tsv'
    prop_defs = {}
    with open(definition, 'r', encoding='utf8') as definitions:
        def_data = csv.reader(definitions, delimiter="\t")
        for row in def_data:
            prop_defs[row[0].rstrip()] = row[4]
    if options.multi is True:
        argos_schemas = {}
        with open(options.input, 'r', encoding='utf8') as file:
            data = csv.reader(file, delimiter="\t")
            header = next(data)
            for row in data:
                if row[1] not in argos_schemas:
                    argos_schemas[row[1]] = {
                        'definitions': {},
                        '$schema': 'http://json-schema.org/draft-07/schema#',
                        '$id': 'https://data.argos.org/schema/'+row[1],
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
                    'examples': row[7],
                    'pattern': row[8]
                }
                if row[2] == 'required':
                    (argos_schemas[row[1]]['required']).append(row[0])
        jsonf = json.dumps(argos_schemas, indent=4)

    else:
        argos_schema = {
            'definitions': {},
            '$schema': 'http://json-schema.org/draft-07/schema#',
            '$id': 'https://data.argos.org/schema/'+ options.input, 
            'title': options.input,
            'type': 'object',
            'required':[],
            'properties':{}
        }

        with open(options.input, 'r', encoding='utf8') as file:
            data = csv.reader(file, delimiter="\t")
            header = next(data)
            for row in data:
                argos_schema['properties'][row[0]] = {
                    '$id':row[3],
                    'title': row[4],
                    'description': row[9],
                    'type': row[5],
                    'default': row[6],
                    'examples': row[7],
                    'pattern': row[8]
                }
                if row[2] == 'required':
                    argos_schema['required'].append(row[0])

        jsonf = json.dumps(argos_schema, indent=4)

    if options.output:
        print(options.output)
        file_name = options.output
        with open(file_name, 'w', encoding='utf-8') as file:
            file.write('[')
            file.write(jsonf)
            file.write(']')
    else:
        print(jsonf)

def main():
    """
    Main function
    """

    usr_args()

if __name__ == "__main__":
    main()
