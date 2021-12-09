import json, csv
import argparse
import jsonschema
from jsonschema import validate
import sys
__version__ = "1.1.0"
__status__ = "Production"

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

    parser.add_argument('-c', '--csv',
                                help="CSV file to create.")

    parser.add_argument('-s', '--schema',
                                required=True,
                                # type = argparse.FileType('r'),
                                help="Root json schema to parse")

    # Print usage message if no args are supplied.
    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    options = parser.parse_args()
    return options

def get_schema(options):
    """
        This function loads the given schema available
    """

    schema_data = {}
    with open(options.schema, 'r') as file:
        schema = json.load(file)
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

def make_csv(options, schema_data): 
    """
        make csv
    """
    if options.csv: 
        csv_file = options.csv
    else: 
        csv_file = 'test.csv'
    csv_columns = ['property', '$id', 'title', 'description', 'type', 'default', 'examples', 'pattern', 'required']
    with open(csv_file, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(csv_columns)
        for key in schema_data:
            row = [key]
            for item in schema_data[key]:
                row.append(schema_data[key][item])
            writer.writerow(row)


def main():
    """
    Main function
    """

    options = usr_args()
    schema_data = get_schema(options)
    make_csv(options, schema_data)

#______________________________________________________________________________#
if __name__ == "__main__":
    main()