#!/usr/bin/env python3
"""ARGOS Data Sheet Validation

This script will perform validation validaton operations for ARGOS data sheets.
General help below.

usage: data_sheet_validator -i INPUT -s SCHEMA [-o OUTPUT] [-v] [-h]

Data Sheet Validation. Used to test a data sheet against a JSON schema. If no
schema is supplied will throw an error

required arguments:
  -i INPUT, --input INPUT
                        Data sheet to validate. This can be a 'tsv' or 'csv'.
                        Other file types are not permitted.
  -s SCHEMA, --schema SCHEMA
                        Root json schema to validate against.

optional arguments:
  -o OUTPUT, --output OUTPUT
                        Output file to create. Default is a JSON file.
  -v, --version         show program's version number and exit
  -h, --help            show this help message and exit
"""

import csv
import json
from argparse import ArgumentParser, SUPPRESS
import os
import sys
from urllib.parse import urlparse
import urllib
import jsonref
from jsonschema import Draft7Validator

__version__ = "0.7"
__status__ = "Draft"

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
        description="Data Sheet Validation. "
            "Used to test a data sheet against a JSON schema. "
            "If no schema is supplied will throw an error")

    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-i', '--input',
        required=True,
        help="Data shet to validate. This can be a 'tsv' or 'csv'. "
            "Other file types are not permitted."
            "If the 'multi' flag is supplied then this should be a directory"
            "path for a group of files" )
    required.add_argument('-s', '--schema',
        required=True,
        help="Root json schema to validate against."
        "If the 'multi' flag is supplied then this should be a directory"
        "path for a group of schema files")

    optional.add_argument('-o', '--output',
        help="Output file to create. Default is a JSON file. ")
    optional.add_argument('-m', '--multi',
        action='store_true',
        help='Flag to indicate if multiple items are being processed')
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

def sheet_2_json(file_path):
    """Convert Data to JSON

    Parameters
    ----------
    file_path: str
        An inpit file to convert to JSON. This should be a TSV/CSV.

    Returns
    -------
        List of JSON objects, one representing each row of the data file
    """

    sheet = []

    extension = file_path.split('.')[-1]
    if extension == 'tsv':
        delimiter='\t'
    elif extension == 'csv':
        delimiter=','
    else:
        print("Must supply a valid 'TSV' or 'CSV'",\
            "\nPlease check your file and try again.")
        exit()

    if os.path.exists(file_path):
        print("Local file supplied")
        with open(file_path, 'r', encoding='utf-8-sig') as file:
            data = csv.reader(file, delimiter=delimiter)
            header = next(data)
            for row in data:
                sheet.append(row)

    elif url_valid(file_path) is True:
        print("Remote file supplied")
        response = urllib.request.urlopen(file_path)
        lines = [l.decode('utf-8') for l in response.readlines()]
        data = csv.reader(lines, delimiter=delimiter)
        header = next(data)
        for row in data:
            sheet.append(row)

    else:
        print("Could not load file. Exiting")
        return

    data_list = []
    for row in sheet:
        line = {}
        for count, item in enumerate(header):
            line[item] = row[count]
        data_list.append(line)
    json_list = json.dumps(data_list)
    return json_list

def validate_schema(in_file, schema_file, out_file):
    """Checks for Schema Complience

    Input and schema are required, output file is optional. If no output is
    supplied then the results will print to the terminal.

    Parameters
    ----------
    in_file: str
        data sheet to be validated

    schema_file: str
        Root json schema to validate against.

    out_file: str, optional
        output file to write errors to. If not provided errors are output to
        the terminal.
    """

    count = 0
    error_flags = 0
    error_strings = {}
    error_strings['input'] = in_file
    error_strings['schema'] = schema_file
    json_list = json.loads(sheet_2_json(in_file))
    no_lines = len(json_list)
    print(f'File with {no_lines} lines supplied.')
    error_strings['lines'] = no_lines

    if os.path.exists(schema_file):
        print("Local schema file supplied")
        with open(schema_file, 'r', encoding='utf-8-sig') as json_schema:
            schema = json.load(json_schema)
    elif url_valid(schema_file) is True:
        print("Remote schema file supplied")
        schema = jsonref.load_uri(schema_file)
    else:
        print("Could not load schema. Exiting")
        return
    error_strings['errors'] = []
    for line in json_list:
        count += 1
        line_errors = {
            'line_number': count,
            'failed_cells': []
        }
        # try:
        validate = Draft7Validator(schema)
        errors = validate.iter_errors(line)

        for item in errors:
            try:
                line_errors['failed_cells'].append(f'{item.relative_path[0]}: {item.message}')
            except IndexError:
                line_errors['failed_cells'].append(f'top_level: {item.message}')
        if len(line_errors['failed_cells']) > 0:
            error_strings['errors'].append(line_errors)
            error_flags += 1

    print(f'{error_flags} lines failed out of {no_lines} lines.')
    if error_flags == 0:
        error_strings['errors'] = 'NONE. Data sheet valid'
    error_strings['summary'] = f'{error_flags} lines failed out of {no_lines}.'

    if out_file:
        with open(out_file, 'w', encoding='utf-8') as file:
            json.dump(error_strings, file, sort_keys = True, indent = 4,
               ensure_ascii = False)
    else:
        print('Input: ', error_strings['input'] )
        print('Schema: ', error_strings['schema'])
        for error in error_strings['errors']:
            try:
                print('Line number', error['line_number'])
                for cell in error['failed_cells']:
                    print('\t',cell)
            except TypeError:
                print(error_strings['errors'])
                break

def file_schema_matcher(input_dir, schema_dir, output_dir):
    """Schema and Datasheet matcher

    Invoked with the multi flag.
    """
    if output_dir is None:
        raise Exception("Output directory required for Multi operation")
    input_list = []
    schema_list = []
    for (root, dirs, file) in os.walk(input_dir):
        for f in file:
            input_list.append(f'{input_dir}/{f}')
    core_dir = f'{schema_dir}/core/'
    non_core_dir = f'{schema_dir}/annotation/'
    
    for (root, dirs, file) in os.walk(core_dir):
        for f in file:
            schema_list.append(f'{core_dir}{f}')
    for (root, dirs, file) in os.walk(non_core_dir):
        for f in file:
            schema_list.append(f'{non_core_dir}{f}')
    # for sheet in input_list:
    #     datasheet = ''
    #     schema = ''
    #     result_file = ''
    count = 0
    tracking_list = []
    for item in schema_list:
        if 'uniprot-proteome' in item:
            search_term = 'uniprot-proteome'
        else:
            search_term = item.split('/')[-1].split('.')[0]
        for sheet in input_list:
            if search_term in sheet:
                count += 1
                output_file = output_dir + '/' + sheet.split('/')[-1].split('.')[0]+'.json'
                print('\n',search_term, [sheet, item])
                validate_schema(sheet, item, output_file)
                tracking_list.append(sheet)

    print('\n')
    if count != len(input_list):
        print(f'{count} files validated', '!=', len(input_list), 'files supplied')
        x = set(input_list)
        y = set(tracking_list)
        print("\nUnvalidated files:")
        for i in x.difference(y):
            print('\t',i)

def main():
    """
    Main function
    """

    options = usr_args()
    if options.multi is True:
        print('Multi is True')
        file_schema_matcher(options.input, options.schema, options.output)
    else:
        validate_schema(options.input, options.schema, options.output)

if __name__ == "__main__":
    main()
