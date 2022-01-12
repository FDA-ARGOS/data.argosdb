#!/usr/bin/env python3
""" SRA Sample Parser

This script will parse the XML result file from a list of SRA samples(SRS). It takes a
directory as an input and will iterate through all files, parsing each one to first a dictionary, and then
outputting to the terminal in TSV format, or to a file if an output file is provided. 

"""

import csv, os
import requests
import xml.etree.ElementTree as ET
__version__ = "0.1.0"
__status__ = "Beta"

xmlfile = 'home/biosampleTest.xml'

trace_url = ['https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=', '&retmode=xml']

samples = {}
header = []

def parseXML( file, dict):
    """
    
    """
    
    biosample_id = sra_id = strain = organism_name = sample_name = ncbi_taxonomy_id = isolate = collected_by = collection_date = geo_loc_name = \
        isolation_source = lat_lon = culture_collection = host = host_age = host_description = host_disease = host_disease_outcome = \
        host_disease_stage = host_health_state = host_sex = identification_method = type_material = '-'

    attribute_list = [strain, isolate, collected_by, collection_date, geo_loc_name, isolation_source, lat_lon, culture_collection, host, host_age, \
        host_description, host_disease, host_disease_outcome, host_disease_stage, host_health_state, host_sex, identification_method, type_material]

    samples = dict
    
    # create element tree object and get root element
    root = ET.parse(file).getroot()
    for biosample in root[:1]:
        biosample_id = biosample.attrib['accession']
        for item in biosample.findall('./'):
            for feature in item.findall('./'):
                # print(feature.tag, feature.text, feature.attrib)
                if 'db' in feature.attrib.keys() and feature.attrib['db'] == 'SRA':
                    sra_id = feature.attrib['db']
                if 'db_label' in feature.attrib.keys():
                    sample_name = feature.text
                if feature.tag == 'Organism':
                    ncbi_taxonomy_id = feature.attrib['taxonomy_id']
                    organism_name = feature.attrib['taxonomy_name']
        for attribute in biosample.findall('./Attributes/'):
            print(attribute.attrib['attribute_name'])
            if attribute.attrib['attribute_name'] == 'strain':
               strain = attribute.text
            if attribute.attrib['attribute_name'] == 'isolate':
               isolate = attribute.text
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


        samples[biosample_id] = [sra_id, organism_name, strain, sample_name, ncbi_taxonomy_id, isolate, collected_by, collection_date, geo_loc_name, \
            isolation_source, lat_lon, culture_collection, host, host_age, host_description, host_disease, host_disease_outcome, host_disease_stage, \
                host_health_state, host_sex, identification_method, type_material]

    print(samples)
    return samples

def printToTerm( dict, list):
    """
    """

    print('\t'.join(item for item in list))
    for run in dict:
        print(run, '\t', '\t'.join(str(item) for item in dict[run]))
    
def main():
    parseXML(xmlfile, samples)
    # printToTerm(samples, header)
      
if __name__ == "__main__":
    main()