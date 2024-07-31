# Updated for schema v1.6 on July 24, 2024

'''This code is used to formulate the tsv for the biosample tsv that is pushed to argos. A text file containing the biosample IDs is needed.This should work on both MacOS and PC computers'''
'''The Python libraries required are in biosample_datagrabber_requirements.py in GITHUB, but here they are:
biopython==1.79
certifi==2022.12.7 # might only be needed for Mac?
numpy==1.22.0
xmltodict==0.13.0'''
#
'''Example Usage:
  python biosample_datagrabber.py --email=myemail@gwu.edu --api_key=myapikey123abc --idfile=mybiosampleids.txt --bco_id=ARGOS_000028'''


from Bio import Entrez
import xmltodict
import time
import argparse

sep = '\t'
# Note that this sleeptime is if authentication (a token) is provided. Use 0.34 if not authenticated.
# Suggesting 0.11 instead of 0.1 just to be safe! Getting banned by NCBI is a huge pain.
sleeptime_withtoken = 0.11 # seconds
sleeptime_notoken = 0.34 # seconds
argos_schema_version = 'v1.6'
schema_keys = ['organism_name',
               'lineage',
               'genome_assembly_id',
               'taxonomy_id',
               'bco_id',
               'schema_version',
               'bioproject',
               'biosample',
               'sra_run_id',
               'instrument',
               'strain',
               'sample_name',
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


    #Get the genome assembly id -----------------------------------------------------------------------
    search = Entrez.esearch(db = 'nucleotide', term = bs_term, retmode='xml', idtype="acc")
    time.sleep(sleeptime)

    record = Entrez.read(search)
    time.sleep(sleeptime)
    #print(f"LINE 101: {record}\n")
    
    assembly_id = "-"
    if record["IdList"]:
        assembly_record_id = record['IdList'][0]
        info = Entrez.efetch(db = 'nucleotide', id = assembly_record_id, rettype="gb", retmode="xml")
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        record_dict = record[0]
        for xref in record_dict["GBSeq_xrefs"]:
            if xref.get("GBXref_dbname") == "Assembly":
                assembly_id = xref["GBXref_id"]

 #Assigning key pair values _________________________________________________________________________________________________
    attr_set['biosample'] = bs_term
    #attr_set['infraspecific_name'] = infra_id #could not get it to work
    attr_set['organism_name'] = sd_json['Description']['Organism']['OrganismName']
    attr_set['genome_assembly_id'] = assembly_id
    attr_set['taxonomy_id'] = sd_json['Description']['Organism']['@taxonomy_id']
    attr_set['bco_id'] = bco_id
    attr_set['schema_version'] = argos_schema_version
    attr_set['bioproject'] = sd_json['Links']['Link']['@label']
    attr_set['sra_run_id'] = SRA_id
    attr_set['instrument'] = 'TBD' # able to successfully get this for each SRA
    attr_set['sample_name'] = sd_json['Description']['Title']
    attr_set['id_method'] = idm_id
    #attr_set['id_method'] = attr_set['identification method']

    search = Entrez.esearch(db = 'sra', term = SRA_id, retmode='xml')
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    sra_id_list = record['IdList']

    bs_data_list = []
    for sra_id in sra_id_list:
        bs_data = attr_set.copy()
        info = Entrez.esummary(db='sra',id=sra_id)
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)
        r = record[0]
        run_info = r['Runs']
        srr_num = xmltodict.parse(run_info)['Run']['@acc']
        # Enclose the XML in an outer set of tags; otherwise the XML is a series of
        # XML objects, not a single one
        exp_xml = "<biosample>" + r['ExpXml'] + "</biosample>"
        exp_json = xmltodict.parse(exp_xml)['biosample']
        # extract the value part of the key:value pair
        bs_data['instrument'] = list(exp_json['Instrument'].values())[0] # FIX THIS, "It works now" - Christie
        bs_data['sra_run_id'] = srr_num
        bs_data_list += [bs_data]

    return bs_data_list

def listify(d, key_order):
    l = []
    for key in key_order:
        l += [d.get(key) or '-']
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
    args = parser.parse_args()

    Entrez.email = args.email
    Entrez.api_key = args.api_key

    output_file = open("biosample.tsv", "w+")
    #print("LINE 187:")
    print (sep.join(schema_keys),
           file=output_file)
    #print()
    for bs_term in args.idfile:  
        bs_term = bs_term.rstrip() # removes any newline character at the end
        #print(f"LINE 191: {bs_term}")
        print(bs_term)
        sleeptime = sleeptime_notoken if args.api_key is None else sleeptime_withtoken
        bs_data_list = bsDataGet(bs_term, sleeptime = sleeptime, bco_id = args.bco_id or '-')

        for bs_data in bs_data_list:
            # read a line from the idfile
            # print its biosample info lines(plural)
            #print("LINE 199:")
            print(sep.join(listify(bs_data, key_order=schema_keys)),
                  file=output_file)
            print()
    #print(bs_term, sleeptime, bco_id)
    output_file.close()
    args.idfile.close()
