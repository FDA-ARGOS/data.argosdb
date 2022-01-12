#!/usr/bin/env python3
"""XML feature parser[WIP]

This script will parse an XML
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
        required=True,
        help="Tag or featured for extraction from XML file.")
    parser.add_argument('-o', '--output',
        required=False,
        help="Output file. If no output is provided, the results will output to the terminal.")
    if len(sys.argv) <= 1:
        sys.argv.append('--help')
    
    return parser.parse_args()

def parseXML(xmlFile, tag):
    """Parse XML file
    
    Parameters
    ----------
    xmlFile: str
        file path/name to be parsed
    tag: str
        tag or featured for extraction from XML file
    """

    tag_list = tag.split('/')
    
    try:
        count = 0
        path = './'
        root = ET.parse(xmlFile).getroot()
        # node = root
        # for index in range(len(tag_list)):
        #     while count <= len(tag_list):
        #         print(index, tag_list[index], count)
        #         for item in node.findall(path):
        #             print(item.tag)
        #         count += 1
        for item in root.findall('./DocumentSummarySet'):
            for run in item.findall('./DocumentSummary'):
                for exp in run.findall('./Runs/'):
                    print(exp.attrib['acc'])
                count += 1
            
        # print(count)
    except ET.ParseError:
        print(xmlFile, 'not well-formed')

def main():
    args = usr_args()
    parseXML(args.file, args.tag)
if __name__ == "__main__":
    main()