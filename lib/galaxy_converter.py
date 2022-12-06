#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Galaxy BCO Converter

Used to convert all Galaxy exported BCO BEFORE Galaxy version 22.09

usage: galaxy_converter.py [options]

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -b BCO, --bco BCO     bco file to convert.
  -e EXTENSION, --extension EXTENSION
                        dataset extension file to add.
  -o OUTPUT, --output OUTPUT
                        File to write output to. If none is supplied, output will print to terminal.
Running the Example:
$ python lib/galaxy_converter.py -b data_files/draft_galaxy_BCO_example.json -e data_files/extsion_example.json -o data_f
iles/final_galaxy_bco_result.json
"""

import json
import argparse
import sys
from tempfile import mkstemp
from shutil import move
from os import fdopen, remove
from urllib.parse import urlparse
from dateutil.parser import parse, ParserError

__version__ = "1.1.0"
__status__ = "Production"


def usr_args():
    """
    functional arguments for process
    """

    parser = argparse.ArgumentParser()

    # set usages options
    parser = argparse.ArgumentParser(prog="galaxy_converter.py", usage="%(prog)s [options]")

    # version
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s " + __version__
    )

    parser.add_argument("-b", "--bco", help="bco file to convert.")
    parser.add_argument("-e", "--extension", help="dataset extension file to add.")
    parser.add_argument(
        "-o",
        "--output",
        # required=True,
        # type = argparse.FileType('r'),
        help="File to write output to. If none is "
        "supplied, output will print to terminal.",
    )
    # Print usage message if no args are supplied.
    if len(sys.argv) <= 1:
        sys.argv.append("--help")

    options = parser.parse_args()
    return options


def timezone_aware(date_text):
    """Makes a datetime Timezone aware

    Parameters
    ----------
    date_text: str
        Datetime object as a string.

    Returns
    -------
    timezone_text: str
        Datetime object as a string with timezone info.
    """
    date_obj = parse(date_text)
    if date_obj.tzinfo is None:
        # print('converting')
        date_t = date_obj.strftime("%Y-%m-%dT%H:%M:%S-04:00")
        return date_t
    # print('no conversion needed')
    return date_text


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


def convert_bco(bco, extension=None):
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

    with open(bco, "r", encoding="utf_8") as file:
        bco_dict = json.load(file)

    if extension is not None:
        bco_dict = add_dataset_extension(bco_dict, extension)

    bco_dict["provenance_domain"]["created"] = timezone_aware(
        bco_dict["provenance_domain"]["created"]
    )
    bco_dict["provenance_domain"]["modified"] = timezone_aware(
        bco_dict["provenance_domain"]["modified"]
    )
    for key in bco_dict["error_domain"]:
        if bco_dict["error_domain"][key] == []:
            bco_dict["error_domain"][key] = {}
    if bco_dict["provenance_domain"]["version"] != str():
        bco_dict["provenance_domain"]["version"] = str(
            bco_dict["provenance_domain"]["version"]
        )
    if bco_dict["provenance_domain"]["review"]:
        del bco_dict["provenance_domain"]["review"]
    for item in bco_dict["provenance_domain"]["contributors"]:
        del item["orcid"]

    for index, item in enumerate(bco_dict["execution_domain"]["script"]):
        if isinstance(item, dict) is False:
            bco_dict["execution_domain"]["script"][index] = {"uri": {"uri": item}}
    if "script_access_type" in bco_dict["execution_domain"]:
        del bco_dict["execution_domain"]["script_access_type"]
    for index, item in enumerate(
        bco_dict["execution_domain"]["external_data_endpoints"]
    ):
        if item is None:
            del bco_dict["execution_domain"]["external_data_endpoints"][index]
    for index, item in enumerate(
        bco_dict["execution_domain"]["software_prerequisites"]
    ):
        if item["uri"]["uri"].startswith("toolshed."):
            bco_dict["execution_domain"]["software_prerequisites"][index]["uri"][
                "uri"
            ] = ("https://" + item["uri"]["uri"])
    for index, item in enumerate(bco_dict["parametric_domain"]):
        if item["step"] != str():
            bco_dict["parametric_domain"][index]["step"] = str(item["step"])
        if item["value"] != str():
            bco_dict["parametric_domain"][index]["value"] = str(item["value"])
    for item in bco_dict["execution_domain"]["software_prerequisites"]:
        if not url_valid(item["uri"]["uri"]):
            item["uri"]["uri"] = "https://" + item["uri"]["uri"]
    for step in bco_dict["description_domain"]["pipeline_steps"]:
        for input_obj in step["input_list"]:
            try:
                input_obj["access_time"] = timezone_aware(input_obj["access_time"])
            except ParserError:
                del input_obj["access_time"]

        for output_obj in step["output_list"]:
            try:
                output_obj["access_time"] = timezone_aware(output_obj["access_time"])
            except ParserError:
                del output_obj["access_time"]

        for preq_obj in step["prerequisite"]:
            try:
                preq_obj["access_time"] = timezone_aware(preq_obj["access_time"])
            except ParserError:
                del preq_obj["access_time"]

    for prereq in bco_dict["execution_domain"]["software_prerequisites"]:
        try:
            prereq["uri"]["access_time"] = timezone_aware(prereq["uri"]["access_time"])
        except ParserError:
            del prereq["uri"]["access_time"]

    for index, input_obj in enumerate(bco_dict["io_domain"]["input_subdomain"]):
        if isinstance(input_obj["uri"], str):
            bco_dict["io_domain"]["input_subdomain"][index] = {"uri": input_obj}
            uri_obj = bco_dict["io_domain"]["input_subdomain"][index]
        try:
            uri_obj["uri"]["access_time"] = timezone_aware(
                uri_obj["uri"]["access_time"]
            )
        except ParserError:
            del input_obj["uri"]["access_time"]
        if not isinstance(uri_obj["uri"]["filename"], str):
            bco_dict["io_domain"]["input_subdomain"][index]["uri"]["filename"] = ""
    for output_obj in bco_dict["io_domain"]["output_subdomain"]:
        try:
            output_obj["uri"]["access_time"] = timezone_aware(
                output_obj["uri"]["access_time"]
            )

        except ParserError:
            del output_obj["uri"]["access_time"]

    return bco_dict


def add_dataset_extension(bco_dict, extension):
    """Add Dataset Extension

    This function loads the given bco and then adds a dataset extension from a
    supplied file.

    Parameters
    ----------
    bco_dict: dict
        Dictionary object with a dataset extension to be added.

    extension:
        File path to BioCompute Object JSON.

    Returns
    -------
    bco_dict: dict
        Dictionary object with a dataset extension added.
    """

    with open(extension, "r", encoding="utf_8") as xfile:
        extension = json.load(xfile)

    bco_dict["extension_domain"].insert(0, extension)

    return bco_dict


def substitution(bco):
    """Subtitution

    This function loads the given bco and then using shell commands substitutes
    broken URIs in the BCO for reconstructed ones based on known hardcoded
    inconsistancies between the Galaxy output and IEEE-2791-2020.

    Parameters
    ----------
    bco: str
        File path to BioCompute Object JSON.

    Returns
    -------
    bco: file
        replacement of the original file object with a mended one
    """

    with open(bco, "r", encoding="utf_8") as file:
        bco_dict = json.load(file)
        for extension in bco_dict["extension_domain"]:
            if "galaxy_extension" in extension.keys():
                galaxy_url = extension["galaxy_extension"]["galaxy_url"][:-1]
    pattern = galaxy_url + "history_content"
    replacement = galaxy_url + "/" + "history_content"

    fhand, abs_path = mkstemp()
    with fdopen(fhand, "w", encoding="utf_8") as new_file:
        with open(bco, "r", encoding="utf_8") as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, replacement))
    remove(bco)
    move(abs_path, bco)


def main():
    """Main function"""

    options = usr_args()
    substitution(options.bco)
    if options.extension:
        jsonf = json.dumps(convert_bco(options.bco, options.extension), indent=4)
    else:
        jsonf = json.dumps(convert_bco(options.bco), indent=4)
    if options.output:
        file_name = options.output
        with open(file_name, "w", encoding="utf-8") as file:
            file.write(jsonf)
    else:
        print(jsonf)


if __name__ == "__main__":
    main()
