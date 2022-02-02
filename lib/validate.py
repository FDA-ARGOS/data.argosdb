#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Validate BioCompute Object

Validates a BCO based on the schema supplied by the user. If the user does not
supply a schema the schema in the BCO will be used.
"""
import sys
import json
import argparse
import jsonschema
from jsonschema import validate
__version__ = "1.1.0"
__status__ = "Production"

def usr_args():
    """Program Arguments

    All arguments for process are defined here.
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

    parser.add_argument('-j', '--json',
                               required=True,
                               help="JSON to process.")

    parser.add_argument('-s', '--schema',
                               # type = argparse.FileType('r'),
                               help="Root json schema to validate against.")

    # Print usage message if no args are supplied.
    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    options = parser.parse_args()
    return options

def get_schema(options):
    """Load Schema
    """
    with open(options.schema, 'r', encoding='utf8') as file:
        schema = json.load(file)
        print('loaded schema')
    return schema


def validate_json(options):
    """REF: https://json-schema.org/ """
    # Describe what kind of json you expect.
    execute_api_schema = get_schema(options)
    with open(options.json, 'r', encoding='utf8') as file:
        data = json.load(file)
        print('loaded data')

    try:
        validate(instance=data, schema=execute_api_schema)
    except jsonschema.exceptions.ValidationError as err:
        print(err)
        err = "Given JSON data is InValid"
        return False, err

    message = "Given JSON data is Valid"
    return True, message

def main():
    """
    Main function
    """

    options = usr_args()
    message = validate_json(options)
    print(message[1])


#______________________________________________________________________________#
if __name__ == "__main__":
    main()
