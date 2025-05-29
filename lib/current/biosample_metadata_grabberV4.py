# Updated for schema v1.6 and the new JSON makeup on May 20, 2025
#Christie Woodside
#Uses the folder of ngsQC files instead of a text file

''' Takes in a folder of ngs JSONs that are from May 15, 2025 and newer. The assembly ID is no lonbger listed at the top of the jsons. This fixes that issue  '''

'''Example Usage:
  python3 biosample_metadata_grabberV4.py --email PERSON@email.gwu.edu --api_key BBBBBBBBB --folder /Users/USER/Desktop/ARGOS/may28/all/ --output /Users/USER/Desktop/ARGOS/may28/BS_clostridiumPerfringens.tsv'''


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
### ADD YOUR API KEYYYYY

apiKey=input("Enter your API key: ")
email=input("Enter your email: ")

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

    #This is the section of new code that distinguishes between v3 nd v4
    # — after you do your esummary on biosample and before you build attr_set —
    # fetch any linked assemblies
    link_handle = Entrez.elink(dbfrom="biosample", db="assembly", id=bs_id)
    time.sleep(sleeptime)
    link_data = Entrez.read(link_handle)
    time.sleep(sleeptime)

    # pull the first linked assembly UID (if there is one)
    assembly_uids = []
    if link_data and link_data[0].get("LinkSetDb"):
        assembly_uids = [link["Id"] for link in link_data[0]["LinkSetDb"][0]["Link"]]

    assembly_acc = ""
    if assembly_uids:
        # get the assembly accession
        sum_handle = Entrez.esummary(db="assembly", id=assembly_uids[0])
        time.sleep(sleeptime)
        sum_rec = Entrez.read(sum_handle)
        time.sleep(sleeptime)
        assembly_acc = sum_rec["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]

    # now put it into your attr_set
    #end of new code

    attr = sd_json['Attributes']['Attribute']
    attr_set = {}

    for att in attr:
        if '#text' in att:
            attr_set[att['@attribute_name']] = att['#text']
        else:
            attr_set[att['@attribute_name']] = None

    ids = sd_json['Ids']['Id']
    ids = sd_json['Attributes']['Attribute']
    idm_id = ''
    for id in ids:
        if id.get('@attribute_name') == 'identification method':
            idm_id = id['#text']

    attr_set['biosample'] = bs_term
    attr_set['infraspecific_name'] = ''
    attr_set['organism_name'] = sd_json.get('Description', {}).get('Organism', {}).get('OrganismName', '')
    #attr_set['genome_assembly_id'] = genome_assembly_id
    attr_set["genome_assembly_id"] = assembly_acc
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
    

    attr_set['instrument'] = 'TBD'
    attr_set['id_method'] = idm_id
    attr_set['bco_id'] = ''
    tax = attr_set['taxonomy_id']

    for t in tax:
        stream = Entrez.efetch(db="Taxonomy", id=tax, retmode="xml")
        records = Entrez.read(stream)
        #print(records[0].keys())
        attr_set['lineage'] = records[0]["Lineage"]

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

#    parser.add_argument('--email', help='email address associated with the NCBI account', type=str)
#    parser.add_argument('--api_key', help='API key associated with the NCBI account (optional)', type=str)
    parser.add_argument('--folder', help='folder containing JSON files with biosample information', type=str, required=True)
    parser.add_argument('--bco_id', help='BCO ID (optional)', type=str)
    parser.add_argument('--output', help='Name you want to assign to the output tsv file', type=str)
    args = parser.parse_args()

    Entrez.email = email
    Entrez.api_key = apiKey

    output_file = open(args.output, "w+")
    print(sep.join(schema_keys), file=output_file)

    sleeptime = sleeptime_notoken if apiKey is None else sleeptime_withtoken

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
