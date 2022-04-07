#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Galaxy BCO Converter

Used to convert all Galaxy exported BCO BEFORE Galaxy version 22.05
"""

from encodings import utf_8
import json
import argparse
import sys

__version__ = "1.0.0"
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

    parser.add_argument('-b', '--bco',
                        help="bco file to convert."
    )
    parser.add_argument('-o', '--output',
                        # required=True,
                        # type = argparse.FileType('r'),
                        help="File to write output to. If none is "
                            "supplied, output will print to terminal."
    )
    # Print usage message if no args are supplied.
    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    options = parser.parse_args()
    return options

def convert_bco(bco):
    """Convert BCO

    This function loads the given bco and then converts based on known
    hardcoded inconsistancies between the Galaxy output and IEEE-2791-2020.

    Parameters
    ----------
    bco: str
        File path to BioCompute Object JSON.

    Returns
    -------
    bco_dict: dict
        Converted dictionary object that when output ot JSON will be 
        IEEE-2791-2020 complient
    """


    with open(bco, 'r', encoding=utf_8) as file:
        bco_dict = json.load(file)
    for key in bco_dict['error_domain']:
        if bco_dict['error_domain'][key] == []:
            bco_dict['error_domain'][key] = {}
    if bco_dict['provenance_domain']['version'] != str():
        bco_dict['provenance_domain']['version'] = str(bco_dict['provenance_domain']['version'])
    for item in bco_dict['provenance_domain']['review']:
        item['reviewer']['contribution'] = list([item['reviewer']['contribution']])
        del item['reviewer']['orcid']
    for item in bco_dict['provenance_domain']['contributors']:
        del item['orcid']
    for index, item in enumerate(bco_dict['execution_domain']['script']):
        if isinstance(item, dict) is False:
            bco_dict['execution_domain']['script'][index] = {'uri': {'uri': item}}
    if bco_dict['execution_domain']['script_access_type']:
        del bco_dict['execution_domain']['script_access_type']
    for index, item in enumerate(bco_dict['execution_domain']['external_data_endpoints']):
        if item is None:
            del bco_dict['execution_domain']['external_data_endpoints'][index]
    for index, item in enumerate(bco_dict['parametric_domain']):
        if item['step'] != str():
            bco_dict['parametric_domain'][index]['step'] = str(item['step'])
        if item['value'] != str():
            bco_dict['parametric_domain'][index]['value'] = str(item['value'])

    return bco_dict


def main():
    """Main function
    """

    options = usr_args()
    jsonf = json.dumps(convert_bco(options.bco), indent=4)
    if options.output:
        file_name = options.output
        with open(file_name, 'w', encoding='utf-8') as file:
            file.write(jsonf)
    else:
        print(jsonf)

#______________________________________________________________________________#
if __name__ == "__main__":
    main()
