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

    parser.add_argument('-c', '--csv',
                                required=True,
                                help="CSV file to create.")

    parser.add_argument('-s', '--schema',
                                # type = argparse.FileType('r'),
                                help="Root json schema to parse")

    # Print usage message if no args are supplied.
    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    options = parser.parse_args()
    return options


def get_csv(options): 
    """
        get csv
    """
    schema = {'required':[], 'properties':{}}
    required = []
    with open(options.csv, 'r') as file:
        data = csv.reader(file)
        header = next(data)
        for row in data:
            # schema.append(row)
            schema['properties'][row[0]] = {
                '$id':row[1],
                'title': row[2],
                'description': row[3],
                'type': row[4],
                'default': row[5],
                'examples': row[6],
                'pattern': row[7]
            }
            # if row[8] == 'required':
            #     schema['required'].append(row[0])
            # if row[9] == 'NA':
            #     continue
            # else:
            #     print(row[9])
            #     schema['properties'][row[0]]['properties'] = json.dumps(row[9])

        
        return schema

def make_schema(options, schema):
    """
        This function writes the schema from the given csv
    """
    print(options.schema)
    
    jsonf = json.dumps(schema, indent=4)

    if options.schema:
        file_name = options.schema
        with open(file_name, 'w', encoding='utf-8') as file: 
            file.write(jsonf)
    else:
        #print(jsonf)
        next
    
    with open('schema/site_qc.json', 'r') as jsonf2:
        new_schema = json.load(jsonf2)
        
        for key in new_schema['properties']:
            if key in schema['properties']:
                new_schema['properties'][key]['description'] = schema['properties'][key]['description'] 
                # new_schema['properties'][key]['description'] = schema['properties'][key]['description'] 
                # new_schema['properties'][key]['description'] = schema['properties'][key]['description'] 
                new_schema['properties'][key]['title'] = schema['properties'][key]['title'] 
    
    print(json.dumps(new_schema, indent=4))
    
    return jsonf

def main():
    """
    Main function
    """

    options = usr_args()
    schema = get_csv(options)
    make_schema(options, schema)



#______________________________________________________________________________#
if __name__ == "__main__":
    main()