# June 17, 2024
# Christie Woodside. Acknowledgements: Hadley
# Github handle: cwoodside1278
#
"""json to tsv code and formatted schema v1.6 ngsQC HIVE3. 
This code will take the JSON file QC output for ngsQC from schema v1.6 stored in a local folder and reformat it into a combined tsv. The output tsv will be pasted and aligned with the 
columns/headers for ngsQC_HIVE (BCI ID ARGOS_000019) found in the data.argosdb.org dataset. This code is used for the datapush to ARGOSdb.
- The code takes in a folder of ngsQC json outputs and formats their information and additional information in a table.
"""
'''The command line input is: 
python3 json2tsv-V2_ngsQC.py -t test_ngs.tsv --schema /Users/steve/desktop/argos/june/ngs --email cool@gwu.edu'''


#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json, glob, os
import csv
import argparse
import sys
import time
from Bio import Entrez
import xmltodict
import re
from xml.etree import ElementTree as ET

__version__ = "1.1.0"
__status__ = "Development"

sleeptime_withtoken = 0.11 # seconds
sleeptime_notoken = 0.35 # seconds
argos_schema_version = 'v1.6'
sep = '\t'
# Note that this sleeptime is if authentication (a token) is provided. Use 0.34 if not authenticated.
# Suggesting 0.11 instead of 0.1 just to be safe! Getting banned by NCBI is a huge pain.
#This code was taken from biosample_datagrabber_v2.py

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
    
    parser.add_argument("-t", "--tsv", help="name of tsv file to create with .tsv extension.") # The name of the output tsv file (and include .tsv extension)
    
    parser.add_argument("-s", "--schema", required=True,
        # type = argparse.FileType('r'), 
        help="Root json schema to parse",
    ) # input ngsQC json file directory/folder

    #So NCBI doesn't ban us
    parser.add_argument('--email', help='Email address associated with the NCBI account', type=str)

    #API key to be faster
    parser.add_argument('--api',
                        help='API key associated with the NCBI account (optional)',
                        type=str) 


    # Print usage message if no args are supplied.
    if len(sys.argv) <= 1:
        sys.argv.append("--help")

    options = parser.parse_args()
    Entrez.email = options.email
    Entrez.api_key = options.api

    # Set the API key if provided
    if options.api:
        Entrez.api_key = options.api

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

def listify(d, key_order):
    l = []
    for key in key_order:
        l += [d.get(key) or '-']
    return l

#Methods that are consistent and don't need to be updated are above^^^^^^^^^^
###################################################################################################################################################################################################################

def bsDataGet(bs_term, sleeptime):
    ''' gets Biosample ID'''
    search = Entrez.esearch(db='sra', term=bs_term, retmode='xml')
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    assem = record['IdList']
    
    info = Entrez.esummary(db="SRA", id=assem)
    time.sleep(sleeptime)
    record = Entrez.read(info)
    # Access the ExpXml and Runs fields
    exp_xml = record[0].get('ExpXml', '')
    # Wrap the ExpXml and Runs in a root element for parsing
    wrapped_exp_xml = f"<root>{exp_xml}</root>"

    if wrapped_exp_xml:
        exp_dict = xmltodict.parse(wrapped_exp_xml)
        exp = exp_dict['root']
        #print (exp['Biosample'])
        return exp['Biosample']

def bsMeta(bs_term, sleeptime):
    ''' gets additional biosample information to add to the tsv'''
    search = Entrez.esearch(db = 'biosample', term = bs_term, retmode='xml') #'ESearch searches and retrieves primary IDs (for use in EFetch, ELink and ESummary) and term translations, and optionally retains results for future use in the user’s environment.'
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    bs_id = record['IdList'][0]
    info = Entrez.esummary(db='biosample',id=bs_id)
    time.sleep(sleeptime)

    record = Entrez.read(info)
    time.sleep(sleeptime)
    r = record['DocumentSummarySet']['DocumentSummary'][0]
    sd = r['SampleData']
    sd_json = xmltodict.parse(sd)['BioSample']

    # Get attributes from "Attributes" section
    attr = sd_json['Attributes']['Attribute']
    attr_set = {}
    for att in attr:
        attr_set[att['@attribute_name']] = att['#text']
    
    #Get ID of SRA link
    ids = sd_json['Ids']['Id']
    SRA_id = ''
    for id in ids:
        if id.get('@db') == 'SRA':
            SRA_id = id['#text']

    #Get the ID for indentification method, because some do not have one
    ids = sd_json['Attributes']['Attribute']
    idm_id = ''
    for id in ids:
        if id.get('@attribute_name') == 'identification method':
            idm_id = id['#text']
        else:
            idm_id

    ids = sd_json['Attributes']['Attribute']
    s_id = ''
    for id in ids:
        if id.get('@attribute_name') == 'strain_name_alias': #get strain
            s_id = id['#text']
        else:
            s_id
    
    #Assigning key pair values _________________________________________________________________________________________________
    attr_set['biosample'] = bs_term
    attr_set['organism_name'] = sd_json['Description']['Organism']['OrganismName']
    attr_set['taxonomy_id'] = sd_json['Description']['Organism']['@taxonomy_id']
    attr_set['schema_version'] = argos_schema_version
    attr_set['bioproject'] = sd_json['Links']['Link']['@label']
    attr_set['instrument'] = 'TBD' # able to successfully get this for each SRA
    attr_set['id_method'] = idm_id
    attr_set['strain'] = s_id

    search = Entrez.esearch(db = 'sra', term = SRA_id, retmode='xml')
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    sra_id_list = record['IdList']

    bs_data_list = []
    for x in sra_id_list:
        bs_data = attr_set.copy()
        info = Entrez.esummary(db='sra',id=x)
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        r = record[0]
        # Enclose the XML in an outer set of tags; otherwise the XML is a series of
        # XML objects, not a single one
        exp_xml = "<biosample>" + r['ExpXml'] + "</biosample>"
        exp_json = xmltodict.parse(exp_xml)['biosample']
        # extract the value part of the key:value pair
        bs_data['instrument'] = list(exp_json['Instrument'].values())[0]
        bs_data['strategy'] = list(exp_json['Library_descriptor'].values())[1]
        bs_data_list += [bs_data]

    return bs_data_list

#from biosample_datagrabber_v2.py on GitHub but used to grab assembly genome accession
def getAssembly(as_term, sleeptime): #removed bco_id value
    '''Get the genome assembly id -----------------------------------------------------------------------'''
    #print(f' as_term: {as_term}')
    search = Entrez.esearch(db='nucleotide', term=as_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)

    if record["IdList"]:
        assembly_record_id = record['IdList'][0]
        info = Entrez.efetch(db = 'nucleotide', id = assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        #print(f'This is record dict in getAssembly \n {record_dict}')
        
        #This is to get the genome assembly id
        for xref in record_dict["GBSeq_xrefs"]:
            if xref.get("GBXref_dbname") == "Assembly":
                assembly_id = xref["GBXref_id"]
                #print(f" this is assembly_id in getAssem {assembly_id}")
                return assembly_id


#to get the lineage
def getLin(l_term, sleeptime): #removed bco_id value
    '''Get the lineage information'''
    search = Entrez.esearch(db = 'nucleotide', term = l_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    lin_id = ""
    if record["IdList"]:
        assembly_record_id = record['IdList'][0] #copy and pasted
        info = Entrez.efetch(db = 'nucleotide', id = assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        lin_id = record_dict["GBSeq_taxonomy"]
        return lin_id


def extract_srr_id(file_source):
    """Remove the '_#.fasta' or '_#.fastq' suffix from the file source. and just the SRR id"""
    if file_source:
        #print(f"Original file source: {file_source}")  # Debugging line
        # Updated regex pattern to ensure it matches correctly
        srr_id = re.sub(r'(_\d+)?\.(fasta|fastq)$', '', file_source)
        #print(f"Extracted SRR ID: {srr_id}")  # Debugging line
        return srr_id if srr_id else "-"
    return "-"

def make_tsv(options):
    columns_data = json.load(open("./columns_ngs.json", "r"))
    with open(options.tsv, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(columns_data["columns"])

        json_files = glob.glob(os.path.join(options.schema, '*.json'))
        # print(f"Looking for JSON files in: {options.schema}")
        # print(f"Found JSON files: {json_files}")

        for schema in json_files:
            with open(schema, "r") as jsonfile:
                data = json.load(jsonfile)
                for item in data[columns_data["top_level"]]:
                    flat_item = flatten_json(item)
                    row = []

                    # Extract SRR ID
                    srr_id = extract_srr_id(flat_item.get("assembled_genome_acc", ""))
                    bs_id = bsDataGet(srr_id, sleeptime_withtoken)
                    bs_data_list = bsMeta(bs_id, sleeptime_withtoken)

                    # Iterate through each column
                    for key in columns_data["columns"]:
                        if key == "biosample":
                            row.append(bs_id if bs_id else "-")
                        elif key == "genome_assembly_id":
                            assembly_id = getAssembly(bs_id, sleeptime_withtoken)
                            row.append(assembly_id if assembly_id else "-")
                        elif key == "lineage":
                            l_id = getLin(bs_id, sleeptime_withtoken)
                            row.append(l_id if l_id else "-")
                        elif key == "sra_run_id":
                            row.append(srr_id)
                        elif key == "bco_id":
                            row.append("ARGOS_000019") #manually added
                        elif key == "ngs_read_file_source":
                            row.append("SRA")
                        else:
                            # If the key is in the bs_data_list, append the data to the row
                            if key in bs_data_list[0]:  # Assuming bs_data_list has only one item
                                row.append(bs_data_list[0][key])
                            else:
                                # Default behavior for JSON data
                                if key in columns_data["header_map"]:
                                    key = columns_data["header_map"][key]
                                row.append(flat_item.get(key, ""))
                    
                    writer.writerow(row)
                        

def main():
    """
    Main function
    """
    options = usr_args()
    Entrez.email = options.email
    Entrez.api_key = options.api
    #directory = options.schema
    make_tsv(options)


# ______________________________________________________________________________#
if __name__ == "__main__":
    options = usr_args()
    Entrez.email = options.email
    Entrez.api_key = options.api
    main()
