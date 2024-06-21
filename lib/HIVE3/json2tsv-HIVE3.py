#!/usr/bin/env python
# -*- coding: utf-8 -*-
# June 18, 2024
#Christie Woodside
import json, csv
import pandas as pd
import argparse
import jsonschema
from jsonschema import validate
import sys

'''json to tsv code used for the HIVE3 QC Outputs. Works for both nested and not nested json files'''

__version__ = "1.1.0"
__status__ = "Development"

def usr_args():
    """
    functional arguments for process
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

    parser.add_argument('-t', '--tsv',
                                help="tsv file to create.")

    parser.add_argument('-s', '--schema',
                                required=True,
                                #type = argparse.FileType('r'),
                                help="Root json schema to parse")

    # Print usage message if no args are supplied.
    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    options = parser.parse_args()
    return options

def flatten_json(y):
    ''' flattens the JSON because it is nested. Code from stack overflow and ChatGPT'''
    out = {}

    def flatten(x, name=''):
        if isinstance(x, dict):
            for a in x:
                flatten(x[a], name + a + '_')
        elif isinstance(x, list):
            i = 0
            for a in x:
                flatten(a, name + str(i) + '_')
                i += 1
        else:
            out[name[:-1]] = x

    flatten(y)
    return out

def get_schema(options):
    """
        This function loads the given schema available
    """

    schema_data = {}
    with open(options.schema, "r") as file:
        schema = json.load(file) 
  
    
     # Assuming "ngsqc" is the main array you want to process. Can change the schema name if you would like to based on the json file.
    if 'ngsqc' not in schema:
        raise KeyError("The key 'ngsqc' is not in the schema JSON file.")
    
    # Process each item in the "ngsqc" array
    properties = []
        #change the schema name here. Replace ngsqc
    for item in schema['ngsqc']:
        flat_item = flatten_json(item)
        properties.extend(flat_item.keys())

    
    # Remove duplicates while preserving order
    properties = list(dict.fromkeys(properties))
    return properties
    

    for key in schema['properties']:
        print(key)
        if schema['properties'][key]['type'] == 'array':
            schema_data[key] = {
                '$id':schema['properties'][key]['$id'],
                'title': schema['properties'][key]['title'],
                'description': schema['properties'][key]['description'],
                'type': schema['properties'][key]['type'],
                'default': schema['properties'][key]['default'],
                'examples': schema['properties'][key]['items']['examples'],
                'pattern': schema['properties'][key]['items']['pattern']
            }    

        elif schema['properties'][key]['type'] == 'number' or schema['properties'][key]['type'] == 'integer':
            schema_data[key] = {
                '$id':schema['properties'][key]['$id'],
                'title': schema['properties'][key]['title'],
                'description': schema['properties'][key]['description'],
                'type': schema['properties'][key]['type'],
                'default': schema['properties'][key]['default'],
                'examples': schema['properties'][key]['examples'],
                'pattern': None
            } 

        elif schema['properties'][key]['type'] == 'object': 
            # schema_data[key] = {
            #     '$id':schema['properties'][key]['$id'],
            #     'title': schema['properties'][key]['title'],
            #     'description': schema['properties'][key]['description'],
            #     'type': schema['properties'][key]['type'],
            #     'default': schema['properties'][key]['default'],
            #     'examples': schema['properties'][key]['examples'],
            #     'pattern': schema['properties'][key]['pattern']
            # }
            continue

        else: 
            schema_data[key] = {
                '$id':schema['properties'][key]['$id'],
                'title': schema['properties'][key]['title'],
                'description': schema['properties'][key]['description'],
                'type': schema['properties'][key]['type'],
                'default': schema['properties'][key]['default'],
                'examples': schema['properties'][key]['examples'],
                'pattern': schema['properties'][key]['pattern']
            }    
        
        if key in schema['required']:
            schema_data[key]['required'] = 'required'
        else:
            schema_data[key]['required'] = 'optional'

    return schema_data

def make_tsv(options, schema_data):
    """
        This function writes the data to a tsv file
    """
    with open(options.tsv, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        
        # Write the header row
        writer.writerow(schema_data)

        with open(options.schema, 'r') as jsonfile:
            data = json.load(jsonfile)
            # this allows for the flattened JSON to be created as a tsv
            for item in data['ngsqc']:
                flat_item = flatten_json(item)
                row = []
                for key in schema_data:
                    row.append(flat_item.get(key, ''))
                writer.writerow(row)
        


def main():
    """
    Main function
    """

    options = usr_args()
    schema_data = get_schema(options)
    make_tsv(options, schema_data)

#______________________________________________________________________________#
if __name__ == "__main__":
    main()
