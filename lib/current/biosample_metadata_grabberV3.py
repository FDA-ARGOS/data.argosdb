# Updated for schema v1.6 and the new JSON IDs on Oct 28, 2024
#Christie Woodside
#Uses the folder of ngsQC files instead of a text file


'''This code is used to formulate the tsv for the biosample tsv that is pushed to argos. 
The code takes in a folder of ngsQC files, which now have the Biosample accession listed at the top of the JSON, 
and used that to get the biosample information.'''

'''The Python libraries required are in biosample_datagrabber_requirements.py in GITHUB, but here they are:
biopython==1.79
certifi==2022.12.7 # might only be needed for Mac?
numpy==1.22.0
xmltodict==0.13.0'''
#
'''Example Usage:
  python3 biosample_metadata_grabberV3.py --email christie.woodside@email.gwu.edu --api_key bfbde99c962d228023e8d62a078bdb12d108 --folder /Users/christiewoodside/Desktop/ARGOS/oct_28/all/ --output /Users/christiewoodside/Desktop/ARGOS/oct_28/BS_clostridiumPerfringens.tsv'''


from Bio import Entrez
import xmltodict
import time
import argparse
import os
import json

sep = '\t'
sleeptime = 0.11

sleeptime_withtoken = 0.11  # seconds
sleeptime_notoken = 0.34  # seconds

'''Api key for CW as of Oct 1 bfbde99c962d228023e8d62a078bdb12d108'''

argos_schema_version = 'v1.6'
schema_keys = ['organism_name',
               'infraspecific_name',
               'lineage',
               'taxonomy_id',
               'bco_id',
               'schema_version',
               'bioproject',
               'biosample',
               'strain',
               'sra_run_id',
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




def bsDataGet(bs_term, genome_assembly_id, srr_id, sleeptime, bco_id):
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
    # for att in attr:
    #     #print(att['#text'])
    #     attr_set[att['@attribute_name']] = att['#text']

    for att in attr:
        if '#text' in att:
            attr_set[att['@attribute_name']] = att['#text']
        else:
            attr_set[att['@attribute_name']] = None
    #print("attr_set:", attr_set)

    ids = sd_json['Ids']['Id']
    # SRA_id = ''
    # for id in ids:
    #     if id.get('@db') == 'SRA':
    #         #print("hello")
    #         SRA_id = id['#text']

    ids = sd_json['Attributes']['Attribute']
    #print(ids)
    idm_id = ''
    for id in ids:
        if id.get('@attribute_name') == 'identification method':
            idm_id = id['#text']

    # assembly_id = "-"
    # if record["IdList"]:
    #     assembly_record_id = record['IdList'][0]
    #     info = Entrez.efetch(db='nucleotide', id=assembly_record_id, rettype="gb", retmode="xml")
    #     time.sleep(sleeptime)
    #     record = Entrez.read(info)
    #     time.sleep(sleeptime)
    #     record_dict = record[0]
    #     for xref in record_dict["GBSeq_xrefs"]:
    #         if xref.get("GBXref_dbname") == "Assembly":
    #             assembly_id = xref["GBXref_id"]

    #print(attr_set)
    attr_set['biosample'] = bs_term
    attr_set['infraspecific_name'] = ''
    #attr_set['organism_name'] = sd_json['Description']['Organism']['OrganismName']
    attr_set['organism_name'] = sd_json.get('Description', {}).get('Organism', {}).get('OrganismName', '')
    attr_set['genome_assembly_id'] = genome_assembly_id
    attr_set['taxonomy_id'] = sd_json['Description']['Organism']['@taxonomy_id']
    attr_set['bco_id'] = bco_id
    attr_set['schema_version'] = argos_schema_version
    attr_set['sample_name'] = sd_json['Description']['Title']
    
    #to fix if there is no bioproject:
    attr_set['bioproject'] = ""
    if sd_json.get('Links') and sd_json['Links'].get('Link'):
        if isinstance(sd_json['Links']['Link'], list):
            attr_set['bioproject'] = sd_json['Links']['Link'][0].get('@label', '')
        else:
            attr_set['bioproject'] = sd_json['Links']['Link'].get('@label', '')

    if 'isolation source' in attr_set:
        attr_set['isolation_source'] = attr_set['isolation source']
    if 'collected-by' in attr_set:
        attr_set['collected_by'] = attr_set['collected-by']
    if 'lat-lon' in attr_set:
        attr_set['lat_lon'] = attr_set['lat-lon']
    if 'sample name' in attr_set:
        attr_set['sample_name'] = attr_set['sample name']
    

    #attr_set['bioproject'] = sd_json['Links']['Link']['@label']
    #attr_set['sra_run_id'] = SRA_id
    attr_set['instrument'] = 'TBD'
    attr_set['id_method'] = idm_id
    attr_set['bco_id'] = ''

    # Populate lineage using genome assembly ID
    #lineage = getLin(genome_assembly_id, sleeptime)
    #print("lineage: ", lineage)
    tax = attr_set['taxonomy_id']
    #print("tax: ",tax)

    for t in tax:
        stream = Entrez.efetch(db="Taxonomy", id=tax, retmode="xml")
        records = Entrez.read(stream)
        #print(records[0].keys())
        attr_set['lineage'] = records[0]["Lineage"]

    

    # search = Entrez.esearch(db='sra', term=SRA_id, retmode='xml')
    # time.sleep(sleeptime)
    # record = Entrez.read(search)
    # time.sleep(sleeptime)
    # sra_id_list = record['IdList']
    sra_id_list = [s.strip() for s in srr_id.split(';') if s.strip()]

    bs_data_list = []
    for sra_id in sra_id_list:
        bs_data = attr_set.copy()
        #info = Entrez.esummary(db='sra', id=sra_id)
        search = Entrez.esearch(db='sra', term=sra_id, retmode='xml')
        time.sleep(sleeptime)
        search_record = Entrez.read(search)
        time.sleep(sleeptime)

        if not search_record['IdList']:
            print(f"DEBUG: Could not find SRR ID: {sra_id}")
            continue  # skip to the next SRR
        #print(f"DEBUG: Successfully processed SRR ID: {sra_id}")

        uid = search_record['IdList'][0]
        info = Entrez.esummary(db='sra', id=uid)
        time.sleep(sleeptime)
        record = Entrez.read(info)
        time.sleep(sleeptime)

        r = record[0]
        #run_info = r['Runs']
        #srr_num = xmltodict.parse(run_info)['Run']['@acc']
        exp_xml = "<biosample>" + r['ExpXml'] + "</biosample>"
        exp_json = xmltodict.parse(exp_xml)['biosample']
        bs_data['instrument'] = list(exp_json['Instrument'].values())[0]
        bs_data['sra_run_id'] = sra_id
        
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

    parser.add_argument('--email', help='email address associated with the NCBI account', type=str)
    parser.add_argument('--api_key', help='API key associated with the NCBI account (optional)', type=str)
    parser.add_argument('--folder', help='folder containing JSON files with biosample information', type=str, required=True)
    parser.add_argument('--bco_id', help='BCO ID (optional)', type=str)
    parser.add_argument('--output', help='Name you want to assign to the output tsv file', type=str)
    args = parser.parse_args()

    Entrez.email = args.email
    Entrez.api_key = args.api_key

    output_file = open(args.output, "w+")
    print(sep.join(schema_keys), file=output_file)

    sleeptime = sleeptime_notoken if args.api_key is None else sleeptime_withtoken

    # Iterate through each JSON file in the specified folder
    for filename in os.listdir(args.folder):
        if filename.endswith('.json'):
            file_path = os.path.join(args.folder, filename)
            with open(file_path, 'r') as json_file:
                data = json.load(json_file)
                
                # Extract the biosample ID (assuming it's the first key in the JSON)
                #biosample_id = list(data.keys())[0]
                biosample_id = data.get("biosample")
                genome_assembly_id = data.get("assembly")
                srr_id = data.get("shortReads")
                print(biosample_id, srr_id)

                # Fetch biosample data
                bs_data_list = bsDataGet(biosample_id, genome_assembly_id, srr_id, sleeptime=sleeptime, bco_id=args.bco_id or '')
                # print(f"DEBUG: biosample={biosample_id}, srr_id={srr_id}")
                # print(f"DEBUG: Got {len(bs_data_list)} bs_data entries")

                # Write each biosample data to the output file
                for bs_data in bs_data_list:
                    #print("bs_data keys:", bs_data.keys())
                    print(sep.join(listify(bs_data, key_order=schema_keys)), file=output_file)
                    # row = listify(bs_data, key_order=schema_keys)
                    # print("row:", row)  # <-- NEW
                    # print(sep.join(row), file=output_file)

    output_file.close()