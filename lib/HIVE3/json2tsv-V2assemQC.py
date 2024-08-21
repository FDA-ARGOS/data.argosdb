#Newly Updated json2tsv-assemQC.py file as of August 20, 2024
# Updated for schema v1.6 August 20, 2024
# Christie Woodside. Acknowledgements: Hadley
# Github handle: cwoodside1278
#
"""json to tsv code and formatted schema v1.6 assemblyQC HIVE3. (qcAll.json)
This code will take the JSON file QC output for assemblyQC from schema v1.6 stored in a local folder and reformat it into a combined tsv. The output tsv will be pasted and aligned with the 
columns/headers for assemblyQC_HIVE (BCI ID ARGOS_000012) found in the data.argosdb.org dataset. This code was used for the datapush to ARGOSdb"""
"""only takes in the folder of jsons not one individual json as of June 24, 2024"""



import json, glob, os
import csv
import argparse
import sys
from Bio import Entrez
import re
import time

__version__ = "1.1.0"
__status__ = "Development"

sleeptime_withtoken = 0.11
sleeptime_notoken = 0.34
argos_schema_version = 'v1.6'
sep = '\t'

# Note that this sleeptime is if authentication (a token) is provided. Use 0.34 if not authenticated.
# Suggesting 0.11 instead of 0.1 just to be safe! Getting banned by NCBI is a huge pain.
#This code was taken from biosample_datagrabber_v2.py

'''Example entry: 
python3 json2tsv-assemQC.py -t test_assembly.tsv --schema /Users/christiewoodside/desktop/argos/data_push_june/all --email christie.woodside@gwu.edu'''

def usr_args():
    """
    Functional arguments for process
    """
    parser = argparse.ArgumentParser(prog="argosdb", usage="%(prog)s [options]")
    # schema version. ex: v1.6
    parser.add_argument("-v", "--version", action="version", version="%(prog)s " + __version__)
    # The name of the output tsv file (and include .tsv extension)
    parser.add_argument("-t", "--tsv", help="TSV file to create.")
    # input assemblyQC (qcAll) json file directory
    parser.add_argument("-s", "--schema", required=True, help="Root JSON schema to parse. (Input JSON file)")
    #So NCBI doesn't ban us
    parser.add_argument('--email', help='Email address associated with the NCBI account', type=str)

    if len(sys.argv) <= 1:
        sys.argv.append("--help")

    options = parser.parse_args()
    Entrez.email = options.email
    return options

def flatten_json(y):
    '''Flattening the JSON input'''
    out = {}

    def flatten(x, name=""):
        if isinstance(x, dict):
            for a in x:
                flatten(x[a], name + a + "_")
        elif isinstance(x, list):
            for i, a in enumerate(x):
                flatten(a, name + str(i) + "_")
        else:
            out[name[:-1]] = x

    flatten(y)
    return out

def bsDataGet(as_term, sleeptime):
    search = Entrez.esearch(db='nucleotide', term=as_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    
    if record["IdList"]:
        assembly_record_id = record['IdList'][0]
        info = Entrez.efetch(db='nucleotide', id=assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        
        for xref in record_dict["GBSeq_xrefs"]:
            if xref.get("GBXref_dbname") == "Assembly":
                return xref["GBXref_id"]
    return ""

def getOrg(o_term, sleeptime):
    search = Entrez.esearch(db='nucleotide', term=o_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    
    if record["IdList"]:
        assembly_record_id = record['IdList'][0]
        info = Entrez.efetch(db='nucleotide', id=assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        return record_dict.get("GBSeq_organism", "")
    return ""

def getTax(t_term, sleeptime):
    search = Entrez.esearch(db='nucleotide', term=t_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    
    if record["IdList"]:
        assembly_record_id = record['IdList'][0]
        info = Entrez.efetch(db='nucleotide', id=assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        
        for feature in record_dict["GBSeq_feature-table"]:
            for qualifier in feature['GBFeature_quals']:
                if qualifier['GBQualifier_name'] == 'db_xref':
                    return qualifier['GBQualifier_value'].replace("taxon:", "")
    return ""

def getLin(l_term, sleeptime):
    search = Entrez.esearch(db='nucleotide', term=l_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    
    if record["IdList"]:
        assembly_record_id = record['IdList'][0]
        info = Entrez.efetch(db='nucleotide', id=assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        return record_dict.get("GBSeq_taxonomy", "")
    return ""

def getGene(g_term, sleeptime):
    """
    Fetches the total number of genes from the genome annotation data.
    """
    search = Entrez.esearch(db='nucleotide', term=g_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)

    if record["IdList"]:
        assembly_record_id = record['IdList'][0]
        info = Entrez.efetch(db='nucleotide', id=assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        
        comment = record_dict.get("GBSeq_comment", "")
        if not comment:
            print(f"No comment field found for accession {assembly_record_id}")
            return "None Found"

        annotation_data = re.search(r'##Genome-Annotation-Data-START##(.*?)##Genome-Annotation-Data-END##', comment, re.DOTALL)
        if annotation_data:
            annotation_data_str = annotation_data.group(1).strip()
            parsed_data = {}
            for line in annotation_data_str.split(';'):
                if '::' in line:
                    key, value = line.split('::', 1)
                    parsed_data[key.strip()] = value.strip()
            return parsed_data.get('Genes (total)', 'Not Found')
        else:
            print(f"No annotation data found for accession {assembly_record_id}")
            return "None Found"
    else:
        print(f"No records found for term {g_term}")
        return "None Found"

def make_tsv(options):
    """
    This function writes the data to a tsv file.
    """
    columns_data = json.load(open("./columns_assembly.json", "r"))
    
    def get_data_from_flat_item(flat_item, key):
        """Retrieve the data based on the key."""
        if key == "genome_assembly_id":
            return bsDataGet(flat_item.get("assembled_genome_acc", ""), sleeptime_withtoken) or "-"
        elif key == "organism_name":
            return getOrg(flat_item.get("assembled_genome_acc", ""), sleeptime_withtoken) or "-"
        elif key == "taxonomy_id":
            return getTax(flat_item.get("assembled_genome_acc", ""), sleeptime_withtoken) or "-"
        elif key == "lineage":
            return getLin(flat_item.get("assembled_genome_acc", ""), sleeptime_withtoken) or "-"
        elif key == "num_genes":
            return getGene(flat_item.get("assembled_genome_acc", ""), sleeptime_withtoken) or "-"
        elif key == "schema_version":
            return 'v1.6'
        else:
            return flat_item.get(columns_data["header_map"].get(key, key), "")

    with open(options.tsv, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")

        # Write the header row
        writer.writerow(columns_data["columns"])

        # Process each JSON file
        for schema in glob.glob(os.path.join(options.schema, '*.json')):
            with open(schema, "r") as jsonfile:
                data = json.load(jsonfile)
                for item in data[columns_data["top_level"]]:
                    flat_item = flatten_json(item)
                    row = [get_data_from_flat_item(flat_item, key) for key in columns_data["columns"]]
                    writer.writerow(row)

def main():
    """
    Main function
    """

    options = usr_args()
    Entrez.email = options.email
    make_tsv(options)


# ______________________________________________________________________________#
if __name__ == "__main__":
    options = usr_args()
    Entrez.email = options.email
    main()