# June 17, 2024
# Christie Woodside. Acknowledgements: Hadley
# Github handle: cwoodside1278
#
"""json to tsv code and formatted schema v1.6 ngsQC HIVE3. 
This code will take the JSON file QC output for ngsQC from schema v1.6 stored in a local folder and reformat it into a combined tsv. The output tsv will be pasted and aligned with the 
columns/headers for ngsQC_HIVE (BCI ID ARGOS_000019) found in the data.argosdb.org dataset. This code wis used for the datapush to ARGOSdb"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json, glob, os
import csv
import argparse
import sys

__version__ = "1.1.0"
__status__ = "Development"


def usr_args():
    """
    functional arguments for process
    """
    parser = argparse.ArgumentParser()
    # set usages options
    parser = argparse.ArgumentParser(prog="argosdb", usage="%(prog)s [options]")
    # schema version. ex: v1.6
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + __version__
    )
    # The name of the output tsv file (and include .tsv extension)
    parser.add_argument("-t", "--tsv", help="tsv file to create.")
    # input ngsQC json file directory
    parser.add_argument(
        "-s",
        "--schema",
        required=True,
        # type = argparse.FileType('r'),
        help="Root json schema to parse",
    )

    # Print usage message if no args are supplied.
    if len(sys.argv) <= 1:
        sys.argv.append("--help")

    options = parser.parse_args()
    return options


def flatten_json(y):
    out = {}

    def flatten(x, name=""):
        if isinstance(x, dict):
            for a in x:
                flatten(x[a], name + a + "_")
        elif isinstance(x, list):
            i = 0
            for a in x:
                flatten(a, name + str(i) + "_")
                i += 1
        else:
            out[name[:-1]] = x

    flatten(y)
    return out


def make_tsv(options):
    """
    This function writes the data to a tsv file
    """
    columns_data = json.load(open("./columns_ngs.json", "r")) #must have columns_ngs.json file stored locally on your computer to work. It is specific to each GC json output
    with open(options.tsv, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")

        # Write the header row
        writer.writerow(columns_data["columns"])
        #take in a folder of jsons to be combined together as a tsv
        json_files = glob.glob(os.path.join(options.schema, '*.json')) 

        for schema in json_files:
            with open(schema, "r") as jsonfile:
                data = json.load(jsonfile)

                for item in data[columns_data["top_level"]]:
                    flat_item = flatten_json(item)
                    row = []
                    for key in columns_data["columns"]:
                        if key in columns_data["header_map"]:
                            key = columns_data["header_map"][key]
                        row.append(flat_item.get(key, ""))
                    writer.writerow(row)


def main():
    """
    Main function
    """
    # json_files = glob.glob(os.path.join(directory, '*.json'))

    # for jsonfile in json_files:
    #     make_tsv(jsonfile)

    options = usr_args()
    #directory = options.schema
    make_tsv(options)


# ______________________________________________________________________________#
if __name__ == "__main__":
    main()
