# Updated for schema v1.6 on Aug 27, 2024

'''This code is used to formulate the tsv for the biosample tsv that is pushed to argos. A text file containing the biosample IDs is needed.This should work on both MacOS and PC computers'''
'''The Python libraries required are in biosample_datagrabber_requirements.py in GITHUB, but here they are:
biopython==1.79
certifi==2022.12.7 # might only be needed for Mac?
numpy==1.22.0
xmltodict==0.13.0'''
#
'''Example Usage:
  python3 biosample_datagrabber.py --email=myemail@gwu.edu --api_key=myapikey123abc --idfile=mybiosampleids.txt --output test_BS.tsv'''
from Bio import Entrez
import xmltodict
import time
import argparse

sep = '\t'

sleeptime_withtoken = 0.11  # seconds
sleeptime_notoken = 0.34  # seconds
argos_schema_version = 'v1.6'
schema_keys = ['organism_name',
               'infraspecific_name',
               'lineage',
               'taxonomy_id',
               'bco_id',
               'schema_version',
               'bioproject',
               'biosample',
               'sra_run_id',
               'strain',
               'genome_assembly_id',
               'sample_name',
               'instrument',
               'isolate',
               'collected_by',
               'collection_date',
               'geo_loc_name',
               'isolation_source',
               'lat_lon',
               'culture_collection',
               'host',
               'host_age',
               'host_description',
               'host_disease',
               'host_disease_outcome',
               'host_disease_stage',
               'host_health_state',
               'host_sex',
               'id_method']

def bsDataGet(bs_term, sleeptime, bco_id):
    search = Entrez.esearch(db='biosample', term=bs_term, retmode='xml')
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    bs_id = record['IdList'][0]
    info = Entrez.esummary(db='biosample', id=bs_id)
    time.sleep(sleeptime)

    record = Entrez.read(info)
    time.sleep(sleeptime)
    r = record['DocumentSummarySet']['DocumentSummary'][0]
    sd = r['SampleData']
    sd_json = xmltodict.parse(sd)['BioSample']

    attr = sd_json['Attributes']['Attribute']
    attr_set = {}
    for att in attr:
        attr_set[att['@attribute_name']] = att['#text']
    ids = sd_json['Ids']['Id']
    SRA_id = ''
    for id in ids:
        if id.get('@db') == 'SRA':
            SRA_id = id['#text']

    ids = sd_json['Attributes']['Attribute']
    idm_id = ''
    for id in ids:
        if id.get('@attribute_name') == 'identification method':
            idm_id = id['#text']

    search = Entrez.esearch(db='nucleotide', term=bs_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)

    assembly_id = "-"
    if record["IdList"]:
        assembly_record_id = record['IdList'][0]
        info = Entrez.efetch(db='nucleotide', id=assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        for xref in record_dict["GBSeq_xrefs"]:
            if xref.get("GBXref_dbname") == "Assembly":
                assembly_id = xref["GBXref_id"]

    attr_set['biosample'] = bs_term
    attr_set['infraspecific_name'] = ''
    attr_set['organism_name'] = sd_json['Description']['Organism']['OrganismName']
    attr_set['genome_assembly_id'] = assembly_id
    attr_set['taxonomy_id'] = sd_json['Description']['Organism']['@taxonomy_id']
    attr_set['bco_id'] = bco_id
    attr_set['schema_version'] = argos_schema_version
    attr_set['bioproject'] = sd_json['Links']['Link']['@label']
    attr_set['sra_run_id'] = SRA_id
    attr_set['instrument'] = 'TBD'
    attr_set['sample_name'] = sd_json['Description']['Title']
    attr_set['id_method'] = idm_id
    attr_set['bco_id'] = 'ARGOS_000020'

    # Populate lineage using genome assembly ID
    lineage = getLin(assembly_id, sleeptime) #some samples may not have a genome and therefore no lineage. Look it up on NCBI Taxonomy
    attr_set['lineage'] = lineage

    search = Entrez.esearch(db='sra', term=SRA_id, retmode='xml')
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    sra_id_list = record['IdList']

    bs_data_list = []
    for sra_id in sra_id_list:
        bs_data = attr_set.copy()
        info = Entrez.esummary(db='sra', id=sra_id)
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        r = record[0]
        run_info = r['Runs']
        srr_num = xmltodict.parse(run_info)['Run']['@acc']
        exp_xml = "<biosample>" + r['ExpXml'] + "</biosample>"
        exp_json = xmltodict.parse(exp_xml)['biosample']
        bs_data['instrument'] = list(exp_json['Instrument'].values())[0]
        bs_data['sra_run_id'] = srr_num
        
        # Check and update the infraspecific_name column
        isolate_value = bs_data.get('isolate', '')
        if isolate_value and isolate_value != 'Not Applicable':
            bs_data['infraspecific_name'] = isolate_value
        
        bs_data_list += [bs_data]

    return bs_data_list


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


def listify(d, key_order):
    l = []
    for key in key_order:
        value = d.get(key) or ''
        
        if key == 'isolate':
            if '/' in value:
                value = value.split('/')[-1]
        
        if key == 'infraspecific_name':
            isolate_value = d.get('isolate')
            if isolate_value and isolate_value != 'Not Applicable' and isolate_value != '':
                if '/' in isolate_value:
                    value = isolate_value.split('/')[-1]
                else:
                    value = isolate_value
        
        l.append(value)
    return l


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--email',
                        help='email address associated with the NCBI account',
                        type=str)
    parser.add_argument('--api_key',
                        help='API key associated with the NCBI account (optional)',
                        type=str)
    parser.add_argument('--idfile',
                        help='text file containing one biosample id on each line',
                        type=argparse.FileType('r'))
    parser.add_argument('--bco_id',
                        help='BCO ID (optional)',
                        type=str)
    parser.add_argument('--output',
                        help= 'Name you want to assign to the output tsv file',
                        type=str)
    args = parser.parse_args()

    Entrez.email = args.email
    Entrez.api_key = args.api_key

    output_file = open(args.output, "w+")
    print(sep.join(schema_keys), file=output_file)

    for bs_term in args.idfile:  
        bs_term = bs_term.rstrip()  # removes any newline character at the end
        sleeptime = sleeptime_notoken if args.api_key is None else sleeptime_withtoken
        bs_data_list = bsDataGet(bs_term, sleeptime=sleeptime, bco_id=args.bco_id or '')

        for bs_data in bs_data_list:
            print(sep.join(listify(bs_data, key_order=schema_keys)),
                  file=output_file)

    output_file.close()
    args.idfile.close()
