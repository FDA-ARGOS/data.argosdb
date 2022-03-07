#!/usr/bin/env python3
"""ARGOS dictionary utilities

This script will perform various functions related to processing the ARGOS
data dictionary and associated schemas. General help below.

positional arguments:
{functions,validate,tsv2json}
    functions           List of all available functions and options.
    validate            Validation options. Used to test a data sheet against
                        a JSON schema. If no schema is supplied will throw an
                        error.
    tsv2json            Used to convert a TSV into a JSON schema. If no mapping
                        file is provided, performs default conversions.
    validate_columns    Validates columns in a list of files in a directory using
                        provided column headers

optional arguments:
-h, --help            show this help message and exit
-v, --version         show program's version number and exit
"""

import csv
import json
import argparse
import sys
import os
from urllib.parse import urlparse
import jsonref
import jsonschema
from jsonschema import validate, ValidationError

__version__ = "0.2.0"
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
    parent_parser.add_argument('-d', '--directory',
        help="Directory")
    parent_parser.add_argument('-c', '--column_index',
        help="Column Index")
    parent_parser.add_argument('-o', '--output',
        help="Output file to create")
    parent_parser.add_argument('-s', '--schema',
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

    # Create a write_schema subcommand
    parser_tsv2json = subparsers.add_parser('write_schema',
        parents=[parent_parser],
        help="Used to convert a TSV into a JSNO schema."
            " If no mapping file is provided, performs default conversions.")
    parser_tsv2json.set_defaults(func=list_2_schema)

    # Create a validate_columns subcommand
    parser_validate_columns = subparsers.add_parser('validate_columns',
        parents=[parent_parser],
        help=" Validates columns in a list of files in a directory using provided column headers")
    parser_validate_columns.set_defaults(func=validate_columns)
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
            print(f"Function: '{choice}'")
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
        Root json schema to validate against.

    options.output: str, optional

    Returns
    -------
        


    """

    json_list = json.loads(sheet_2_json(options.input))
    count = 0
    print(len(json_list))
    if options.schema is None:
        print("ERROR! No schema was supplied. Exiting")
        return

    if os.path.exists(options.schema):
        print("Local file supplied")
        with open(options.schema, 'r', encoding='utf-8-sig') as json_schema:
            schema = json.load(json_schema)
        print('schema worked')
    elif url_valid(options.schema) is True:
        print("Remote file supplied")
        schema = jsonref.load_uri(options.schema)

    for line in json_list:
        count += 1
        try:
            validate(instance=line, schema=schema)
            print(f'Item number {count} is VALID')
        except ValidationError as err:
            print(f'Line {count} failed. {err.message}')
            # err = "Given JSON data is InValid"
            # return False, err

def list_2_schema(options):
    """Create Schema JSON

    Takes a TSV and writes JSON formatted schema

    Parameters
    ----------
    options.input: str
        An inpit file to create the schema/schemas. This should be a TSV.

    options.multi: bool
        Default is False. If true the input file is treated as a flat version
        of a multiple schemas. These will output to a list of JSON.

    options.directory: str, optional
        An output file. If this is supplied then the function output will be
        written to this directory.

    definition: str
        a list of property definitions.
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
    print(options.directory)
    if os.path.exists(options.directory):
        argos_schemas = {}
        with open(options.input, 'r', encoding='utf8') as file:
            data = csv.reader(file, delimiter="\t")
            next(data)
            count = 0
            for row in data:
                count+= 1
                try:
                    prop_defs[row[0].rstrip()]
                except KeyError:
                    print(f"Error! {row[0]} at row {count} is not defined in {definition}. Exiting.")
                    return
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
                    'examples': [row[7]],
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
            next(data)
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

    if os.path.exists(options.directory):
        for item in argos_schemas:
            output = options.directory+item.split('.')[0]+'.json'
            print(output)
            file_name = options.directory+item.split('.')[0]+'.json'
            with open(file_name, 'w', encoding='utf-8') as file:
                file.write(json.dumps(argos_schemas[item], indent=4))
    else:
        print('jsonf')

def validate_columns(options):
    """Validate Columns

    Parameters
    ----------
    options.input: str
        An inpit file to create the schema/schemas. This should be a TSV.

    Returns
    -------
    """
    columns = []
    missing_keys = {}
    try:
        int(options.column_index)
    except ValueError:
        print('Column index is not a number')
        return
    count = 0
    reader = open(options.input, 'r', encoding='utf8')
    for i in reader:
        columns.append(i.split("\t")[int(options.column_index)])
        count += 1
    with open('test2.json', 'w', encoding='utf8') as outfile:
        json.dump(columns, outfile)
        print(columns, outfile)
    for filename in os.listdir(options.directory):
        delimiter=""
        if filename.endswith(".tsv"):
            delimiter="\t"
        elif filename.endswith(".csv"):
            delimiter = ","
        else:
            continue
        file_path = os.path.join(options.directory, filename)
        if os.path.isfile(file_path):
            data = open(file_path, 'r', encoding='utf8')
            header=data.readline().split(delimiter)
            header = [x.strip("\n").strip("\"").lower() for x in header]
            missing_keys[filename] = [x for x in header if x not in columns]
        with open(options.output, 'w', encoding='utf8') as outfile:
            print(missing_keys, outfile)
            json.dump(missing_keys, outfile)

def sheet_2_json(file_path):
    """Convert Data to JSON

    Parameters
    ----------
    file_path: str
        An inpit file to convert to JSON. This should be a TSV/CSV.

    Returns
    -------
        List of JSONs, one representing each row of the data file
    """
    extension = file_path.split('.')[-1]
    if extension == 'tsv':
        delimiter='\t'
    elif extension == 'csv':
        delimiter=','
    sheet = []
    with open(file_path, 'r', encoding='utf-8-sig') as file:
        data = csv.reader(file, delimiter=delimiter)
        header = next(data)
        for row in data:
            sheet.append(row)
    data_list = []
    for row in sheet[:10]:
        line = {}
        for count, item in enumerate(header):
            cell = f'{item}: {row[count]},'
            line[item] = row[count]
        data_list.append(line)
    json_list = json.dumps(data_list)
    return json_list

def url_valid(url):
    """Validate a URL

    Parameters
    ----------
    url: str
        String as a URL to validate
    Returns
    -------
    bool
        A True or False value for the URL supplied

    """

    try:
        result = urlparse(url)
        return all([result.scheme, result.netloc])
    except ValueError:
        return False

def main():
    """
    Main function
    """

    usr_args()

if __name__ == "__main__":
    main()
