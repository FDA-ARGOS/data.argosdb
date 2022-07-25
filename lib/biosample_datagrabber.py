#
# The Python libraries required are in biosample_datagrabber_requirements.py
#

from Bio import Entrez
import xmltodict
import time
import argparse

sep = '\t'
sleeptime = 0.1 # seconds
argos_schema_version = 'v0.7'
schema_keys = ['organism_name',
               'lineage',
               'genome_assembly_id',
               'taxonomy_id',
               'fdaargos_id',
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

def bsDataGet(bs_term):

    search = Entrez.esearch(db = 'biosample', term = bs_term, retmode='xml')
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

    # Get ID of SRA link
    ids = sd_json['Ids']['Id']
    SRA_id = ''
    for id in ids:
        if id.get('@db') == 'SRA':
            SRA_id = id['#text']

    # Get the genome assembly id
    search = Entrez.esearch(db = 'assembly', term = bs_term, retmode='xml')
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    assembly_record_id = record['IdList'][0]
    info = Entrez.esummary(db = 'assembly',id = assembly_record_id)
    record = Entrez.read(info)
    genome_assembly_id = record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']

    attr_set['biosample'] = bs_term
    attr_set['organism_name'] = sd_json['Description']['Organism']['OrganismName']
    attr_set['genome_assembly_id'] = genome_assembly_id
    attr_set['taxonomy_id'] = sd_json['Description']['Organism']['@taxonomy_id']
    attr_set['fdaargos_id'] = attr_set['strain']
    attr_set['bco_id'] = ''
    attr_set['schema_version'] = 'v0.7'
    attr_set['bioproject'] = sd_json['Links']['Link']['@label']
    attr_set['sra_run_id'] = SRA_id
    attr_set['instrument'] = 'TBD' # able to successfully get this for each SRA
    attr_set['sample_name'] = sd_json['Description']['Title']
    attr_set['id_method'] = attr_set['identification method']
    attr_set['schema_version'] = argos_schema_version

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
        bs_data['instrument'] = list(exp_json['Instrument'].values())[0] # FIX THIS
        bs_data['sra_run_id'] = srr_num
        bs_data_list += [bs_data]

    return bs_data_list

def listify(d, key_order):
    l = []
    for key in key_order:
        l += [d.get(key) or '']
    return l

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--email',
                        help='email address associated with the Entrez API account',
                        type=str)
    parser.add_argument('--idfile',
                        help='text file containing one biosample id on each line',
                        type=argparse.FileType('r'))
    args = parser.parse_args()

    Entrez.email = args.email
    #print(args.idfile)
    output_file = open("biosample.tsv", "w+")
    print (sep.join(schema_keys),
           file=output_file)
    for bs_term in args.idfile:
        bs_term = bs_term.rstrip() # removes any newline character at the end
        print(bs_term)
        bs_data_list = bsDataGet(bs_term)

        for bs_data in bs_data_list:
            # read a line from the idfile
            # print its biosample info lines(plural)
            print(sep.join(listify(bs_data, key_order=schema_keys)),
                  file=output_file)

    output_file.close()
    args.idfile.close()
