#!/usr/bin/env python3
""" SRA Biosample Parser

This script will parse the XML result file from a list of SRA samples(SRS). It takes a
single XML file as input, parsing it first a dictionary, and then
outputting to the terminal in TSV format, or to a file if an output file is provided.
"""

import os
import sys
import csv
import argparse
import xml.etree.ElementTree as ET

__version__ = "0.1.1"
__status__ = "Beta"

PLACEHOLDER_CHAR = ''

def usr_args():
    """Program Arguments

    All arguments for process are defined here.
    """

        # initialize parser
    parser = argparse.ArgumentParser()

    # set usages options
    parser = argparse.ArgumentParser(
        prog='sra_trace_xml.py',
        usage='%(prog)s [options]')

    # version
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s ' + __version__)

    parser.add_argument('-f', '--file',
        required=True,
        help="Input file. The input file should be a collection of SRS XML \
            files.")
    parser.add_argument('-b', '--bco_id',
        help="BCO ID")
    parser.add_argument('-s', '--schema',
        help="Schema version")
    parser.add_argument('-o', '--output',
        help="Output file. If no output is provided, the resulting table will \
            be output to the terminal.")
    if len(sys.argv) <= 1:
        sys.argv.append('--help')

    return parser.parse_args()

def parse_xml( xml_file, samples, bco_id, schema_version):
    """Parse XML file

    Parameters
    ----------
    xml_file: str
        file path/name to be parsed
    samples: dict
        dictionary of samples

    Returns
    -------
    samples: dict
        dictionary of samples
    """

    organism_name = infraspecific_name = lineage = taxonomy_id \
        = bioproject = biosample \
        = strain = genome_assembly_id = isolate = sample_name \
        = instrument = collected_by = collection_date \
        = geo_loc_name = isolation_source = lat_lon = culture_collection = host\
        = host_age = host_description = host_disease = host_disease_outcome \
        = host_disease_stage = host_health_state = host_sex \
        = id_method = biosample_score = PLACEHOLDER_CHAR

    bco_id = 'ARGOS_000010'
    bioproject = 'PRJNA231221'
    schema_version = schema_version or PLACEHOLDER_CHAR

    # create element tree object and get root element
    root = ET.parse(xml_file).getroot()
    for biosample in root:
        biosample_id = biosample.attrib['accession']
        for item in biosample.findall('./'):
            for feature in item.findall('./'):
                if 'db' in feature.attrib.keys() and feature.attrib['db'] == 'SRA':
                    sra_id = feature.text
                if 'db_label' in feature.attrib.keys():
                    sample_name = feature.text
                if feature.tag == 'Organism':
                    taxonomy_id = feature.attrib['taxonomy_id']
                    organism_name = feature.attrib['taxonomy_name']
         
        for attribute in biosample.findall('./Attributes/'):
            if attribute.attrib['attribute_name'] == 'strain':
                strain = attribute.text
            if attribute.attrib['attribute_name'] == 'isolate':
                isolate= attribute.text
            if attribute.attrib['attribute_name'] == 'sample_name':
                sample_name = attribute.text
            if attribute.attrib['attribute_name'] == 'collected_by':
                collected_by = attribute.text
            if attribute.attrib['attribute_name'] == 'collection_date':
                collection_date = attribute.text
            if attribute.attrib['attribute_name'] == 'geo_loc_name':
                geo_loc_name = attribute.text
            if attribute.attrib['attribute_name'] == 'isolation_source':
                isolation_source  = attribute.text
            if attribute.attrib['attribute_name'] == 'lat_lon':
                lat_lon  = attribute.text
            if attribute.attrib['attribute_name'] == 'culture_collection':
                culture_collection = attribute.text
            if attribute.attrib['attribute_name'] == 'host':
                host = attribute.text
            if attribute.attrib['attribute_name'] == 'host_age':
                host_age  = attribute.text
            if attribute.attrib['attribute_name'] == 'host_description':
                host_description = attribute.text
            if attribute.attrib['attribute_name'] == 'host_disease':
                host_disease = attribute.text
            if attribute.attrib['attribute_name'] == 'host_disease_outcome':
                host_disease_outcome = attribute.text
            if attribute.attrib['attribute_name'] == 'host_disease_stage':
                host_disease_stage = attribute.text
            if attribute.attrib['attribute_name'] == 'host_health_state':
                host_health_state = attribute.text
            if attribute.attrib['attribute_name'] == 'host_sex':
                host_sex = attribute.text
            if attribute.attrib['attribute_name'] == 'identification method':
                identification_method = attribute.text
            if attribute.attrib['attribute_name'] == 'type-material':
                type_material = attribute.text
        

        samples[biosample_id] = [organism_name, infraspecific_name, lineage, taxonomy_id, bco_id, \
            schema_version, bioproject, biosample_id, strain, genome_assembly_id, \
            sample_name, instrument, isolate, collected_by, collection_date, \
            geo_loc_name, isolation_source, lat_lon, culture_collection, host,\
            host_age, host_description, host_disease, host_disease_outcome, \
            host_disease_stage, host_health_state, host_sex, \
            identification_method]

        biosample_score = len(samples[biosample_id])-samples[biosample_id].count(PLACEHOLDER_CHAR)
        samples[biosample_id].append(biosample_score)


    return samples

def sample_output( samples, header, output):
    """Sample Output

    If an output file is supplied in the user arguments the samples dictionary
    will be output in TSV to the supplied file. If no output is provided, the
    samples dictionary will print to the terminal in a TSV format.

    Parameters
    ----------
    header: lst of str
        List of column headers for the output
    samples: dict
        dictionary of samples
    outputs: str, optional
        file path/name to output data to
    """

    if output:
        sample_file = os.path.abspath(output)
        with open(sample_file, 'w',  encoding='utf8', newline = '') as file:
            writer = csv.writer(file, header, delimiter='\t')
            writer.writerow(header)
            for key in samples:
                #row = [key] # biosample_id is no longer first in order
                row = []
                for item in range(len(samples[key])):
                    row.append(samples[key][item])
                writer.writerow(row)
        print('Samples written to ', sample_file)
    else:
        print('\t'.join(item for item in header))
        for run in samples:
            print('\t'.join(str(item) for item in samples[run]))
    


def main():
    """Main Function
    """
    samples = {}
    header = ['organism_name', 'infraspecific_name', 'lineage', 'taxonomy_id', 'bco_id', \
        'schema_version', 'bioproject', 'biosample', 'strain', 'genome_assembly_id',\
        'sample_name', 'instrument', 'isolate', 'collected_by',\
        'collection_date', 'geo_loc_name', 'isolation_source', 'lat_lon',\
        'culture_collection', 'host', 'host_age', 'host_description', \
        'host_disease', 'host_disease_outcome', 'host_disease_stage',\
        'host_health_state', 'host_sex', 'id_method',\
        'biosample_score']
    args = usr_args()
    parse_xml(args.file, samples, args.bco_id, schema_version = args.schema)
    sample_output(samples, header, args.output)

if __name__ == "__main__":
    main()
