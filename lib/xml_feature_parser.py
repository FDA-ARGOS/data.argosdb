#!/usr/bin/env python3
"""XML feature parser[WIP]

This script will parse an XML for a single feature and return a list via terminal
or file, if an output is supplied.
"""
__version__ = "0.2.0"
__status__ = "Production"

import sys
import argparse
import xml.etree.ElementTree as ET

def usr_args():
    """Program Arguments

    All arguments for process are defined here.
    """

        # initialize parser
    parser = argparse.ArgumentParser()

    # set usages options
    parser = argparse.ArgumentParser(
        prog='XML feature parser',
        usage='%(prog)s [options]')

    # version
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s ' + __version__)

    parser.add_argument('-f', '--file',
        required=True,
        help="Path for input file. The input file should be an XML file.")
    parser.add_argument('-t', '--tag',
        help="Tag or featured for extraction from XML file.")
    parser.add_argument('-o', '--output',
        required=False,
        help="Output file. If no output is provided, the results will output to the terminal.")
    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    return parser.parse_args()

def parse_xml(xml_file, output):
    """Parse XML file

    Parameters
    ----------
    xml_file: str
        file path/name to be parsed
    output: str, optional
        file path/name for data to be output to
    """

    # tag_list = tag.split('/')
    items = []
    print('start')
    try:
        count = 0
        root = ET.parse(xml_file).getroot()
        for item in root.findall('.'):
            print("FOUND ONE")
            for run in item.findall('./DocumentSummary'):
                for exp in run.findall('./Runs/'):
                    items.append(exp.attrib['acc'])
                count += 1

        if output:
            with open(output, 'w', encoding='utf-8') as file:
                for item in items:
                    file.write(item+'\n')
        else:
            for item in items:
                print(item)

    except ET.ParseError:
        print(xml_file, 'not well-formed')

def main():
    """Main"""
    args = usr_args()
    parse_xml(args.file, args.output)

if __name__ == "__main__":
    main()
