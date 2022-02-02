#!/usr/bin/env python3
"""XML feature parser[WIP]

This script will parse an XML for a single feature and return a list via terminal
or file, if an output is supplied.
"""
__version__ = "0.1.0"
__status__ = "Beta"
import os
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

def parseXML(xmlFile, output):
    """Parse XML file
    
    Parameters
    ----------
    xmlFile: str
        file path/name to be parsed
    output: str, optional
        file path/name for data to be output to
    """

    # tag_list = tag.split('/')
    items = []
    try:
        count = 0
        root = ET.parse(xmlFile).getroot()
        for item in root.findall('./DocumentSummarySet'):
            for run in item.findall('./DocumentSummary'):
                for exp in run.findall('./Runs/'):
                    items.append(exp.attrib['acc'])
                count += 1
            
        
        if output:
            with open(output, 'w') as file:
                for item in items:
                    file.write(item+'\n')
        else:
            for item in items:
                print(item)

    except ET.ParseError:
        print(xmlFile, 'not well-formed')

def main():
    args = usr_args()
    parseXML(args.file, args.output)
if __name__ == "__main__":
    main()