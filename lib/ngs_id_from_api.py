
import os
import requests
import json
import csv
import pandas as pd
pd.options.mode.chained_assignment = None  # removes warning message from overwriting

api_url = "https://beta-api.argosdb.org/records/search"
HIVE_results = []
Crandall_results = []
HIVE_data = [{
  "bcoid": "ARGOS_000019",
  "offset": 1,
  "limit": 10000
}]
Crandall_data = [{
  "bcoid": "ARGOS_000025",
  "offset": 1,
  "limit": 10000
}]

### Need BCO for Pond NGS Data

# defines function to get HIVE Lab specific data
def get_HIVE():
    for item in HIVE_data:
        response = requests.post(api_url, json=item)
        HIVE_results.append(response.json())
    with open('hiveapi.json', 'w', newline = '', encoding='utf-8') as f:
            json.dump(HIVE_results,f, ensure_ascii=False, indent=4)

#Convert json to tsv and make keys into column headers
    with open('hiveapi.json') as jsonfile:
        data = json.load(jsonfile)
    records=data[0]['recordlist']
    datafile = open('hiveTEMP.tsv','w', newline = '')
    tsvwriter=csv.writer(datafile, delimiter= '\t')
    count=0
    for record in records:
        if count==0:
            header = record.keys()
            tsvwriter.writerow(header)
            count+=1
        tsvwriter.writerow(record.values())
#Read/write tsv - writes appropriate headers for ngs_id then writes corresponding columns from ngsQC file
    with open('hiveTEMP.tsv', 'r') as source:
        reader = csv.reader(source, delimiter='\t')
        with open('hiveapi.tsv', 'w', newline = '') as result:
            writer=csv.writer(result, delimiter='\t')
            headings=next(reader)
            #append headers for ngs_id
            writer.writerow(['organism_name', 'intraspecific_name', 'lineage', 'genome_assembly_id', 'taxonomy_id', 'bioproject','biosample','sra_run_id','ngs_read_file_source','ref_org','isolate_identifiers','selection_notes','lab_name','files_processed'])
            # These columns need to be updated/changed to reflect new columns added to ngs_id
            for r in reader:
                writer.writerow([r[3],r[4],r[5],r[6],r[7],r[13],r[14], r[15],r[17]])
    datafile.close()
get_HIVE()

data_hl = pd.read_table('hiveapi.tsv', sep='\t', on_bad_lines='skip')

df_hl = data_hl
#This line originally filtered assembly ids
df_hl=df_hl[df_hl['sra_run_id'].isnull() | ~df_hl[df_hl['sra_run_id'].notnull()].duplicated(subset='sra_run_id',keep='first')]
df_hl.lab_name = 'HIVE Lab'
df_hl.files_processed = 'ngsQC_HL'

    #for record in response.json()['recordlist']:
    #    if record['genome_assembly_id'] not in assemblies:
    #        assemblies.append(record['genome_assembly_id'])
    #        os.system(f"efetch -db assembly -id {record['genome_assembly_id']} -format docsum > test5/{record['genome_assembly_id']}.xml")
    #print(response.status_code)

def get_Crandall():
    for item in Crandall_data:
        response = requests.post(api_url, json=item)
        Crandall_results.append(response.json())
    with open('Crandallapi.json', 'w', newline = '', encoding='utf-8') as f:
            json.dump(Crandall_results,f, ensure_ascii=False, indent=4)
#Convert json to tsv and make keys into column headers
    with open('Crandallapi.json') as jsonfile:
        data = json.load(jsonfile)
    records=data[0]['recordlist']
    datafile = open('CrandallTEMP.tsv','w', newline = '')
    tsvwriter=csv.writer(datafile, delimiter= '\t')
    count=0
    for record in records:
        if count==0:
            header = record.keys()
            tsvwriter.writerow(header)
            count+=1
        tsvwriter.writerow(record.values())


#Read/write tsv - writes appropriate headers for ngs_id then writes corresponding columns from ngsQC file
#
#
# Need to figure out how to pull assembly id for corresponding Biosample
    with open('CrandallTEMP.tsv', 'r') as source:
        reader = csv.reader(source, delimiter='\t')
        with open('Crandallapi.tsv', 'w', newline = '') as result:
            writer=csv.writer(result, delimiter='\t')
            headings=next(reader)
            #append headers for ngs_id
            ####
            ####  ngsQC_Crandall hasn't been updated for v9 - this will need to be edited
            #####
            writer.writerow(['organism_name', 'intraspecific_name', 'lineage', 'genome_assembly_id', 'taxonomy_id', 'bioproject','biosample','sra_run_id','ngs_read_file_source','ref_org','isolate_identifiers','selection_notes','lab_name','files_processed'])
            # These columns need to be updated/changed to reflect new columns added to ngs_id
            for r in reader:
                writer.writerow([r[3],r[4],r[5],r[6],r[7],r[13],r[14], r[15],r[17]])

    datafile.close()
get_Crandall()


data_c = pd.read_table('Crandallapi.tsv', sep='\t', on_bad_lines='skip')
df_c = data_c
#This line originally filtered assembly ids
df_c=df_c[df_c['sra_run_id'].isnull() | ~df_c[df_c['sra_run_id'].notnull()].duplicated(subset='sra_run_id',keep='first')]
df_c.lab_name = 'Crandall Lab'
df_c.files_processed = 'ngsQC_Crandall'

#def get_Pond():
#    for item in Pond_data:
#        response = requests.post(api_url, json=item)
#        Pond_results.append(response.json())
#    with open('Pondapi.json', 'w', newline = '', encoding='utf-8') as f:
#            json.dump(Pond_results,f, ensure_ascii=False, indent=4)
#
#Convert json to tsv and make keys into column headers
#    with open('Pondapi.json') as jsonfile:
#        data = json.load(jsonfile)
#    records=data[0]['recordlist']
#    datafile = open('Pondapi.tsv','w', newline = '')
#    tsvwriter=csv.writer(datafile, delimiter= '\t')
#    count=0
#    for record in records:
#        if count==0:
#            header = record.keys()
#            tsvwriter.writerow(header)
#            count+=1
#        tsvwriter.writerow(record.values())


#Read/write tsv - writes appropriate headers for ngs_id then writes corresponding columns from ngsQC file
#
#
# Need to figure out how to pull assembly id for corresponding Biosample
#    with open('Pondapi.tsv', 'r+') as source:
#        reader = csv.reader(source, delimiter='\t')

#        with open('Pondapi.tsv', 'r+', newline = '') as result:
#            writer=csv.writer(result, delimiter='\t')
#            headings=next(reader)
            #append headers for ngs_id
#            writer.writerow(['organism_name', 'leaf_node', 'genome_assembly_id', 'taxonomy_id', 'bioproject','biosample','sra_run_id','ngs_read_file_source','ref_org','isolate_identifiers','selection_notes','lab_name','files_processed'])


            # These columns need to be updated/changed to reflect new columns added to ngs_id
#            for r in reader:
#                writer.writerow([r[3],r[4],r[6],r[5],r[11],r[12],r[23],r[13]])



#    datafile.close()

#get_Pond()











dfs = [df_hl, df_c]
#
#### dfs = [df_hl, df_c, df_p]
#
combined_file = pd.concat(dfs)
dfinal=combined_file
#Populate Reference Orgs
def ref_org(row):
    if row['genome_assembly_id'] == 'GCA_000865725.1': # (A/Puerto Rico/8/1934(H1N1))
        return "Yes"
    elif row['genome_assembly_id'] == 'GCA_009858895.3': #isolate Wuhan-Hu-1
        return "Yes"
    elif row['genome_assembly_id'] == 'GCA_001558355.2': #LT2
        return "Yes"
    elif row['genome_assembly_id'] == 'GCA_000857325.2': #Marburg
        return "Yes"
    elif row['genome_assembly_id'] == 'GCA_003102975.1': #HXB2
        return "Yes"
    else:
        return "No"
dfinal=dfinal.assign(ref_org=dfinal.apply(ref_org, axis =1))

#Populate selection_notes for all organisms in ARGOS Bioproject
def selection_notes(row):
    if row['bioproject'] == 'PRJNA231221':
        return "Belongs to FDA-ARGOS PRJNA231221."
    elif row['bioproject'] == 'PRJNA726840':
        return "Coding-complete Genome sequences for SARS-CoV-2 B.1.1.7 and B.1.351 Variants from Metro Manila, Philippines, outlined in the following paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8281087/"
    elif row['bioproject'] ==  'PRJNA729484':
        return "Raw sequencing reads were collected for the SARS-CoV-2 P.1 variant in Northeast Brazil, outlined in this paper:  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8321350/"
    elif row['bioproject'] ==  'PRJNA791622' and row['lab_name'] == 'Pond Lab':
        return "For omicron, we are selecting EPI_ISL_6913953. Sequencing was conducted on  Illumina MiSeq, has high coverage, and a consistent quality score across all base calls above 30. Raw reads are available at https://www.ncbi.nlm.nih.gov/sra/SRX13486794, and a full description of the patient harboring the virus is supplied with the following publication https://academic.oup.com/cid/advance-article/doi/10.1093/cid/ciab1072/6494531?login=true. The patient was one of the first two known    COVID-19 cases classified as omicron in Japan. To put the collection date of 28 November 2021 in perspective, the first known omicron sample was collected on 8 November 2021. Raw reads from South Africa are available, but the average phred quality score is much lower for those samples."
    elif row['bioproject'] ==  'PRJNA791622' and row['lab_name'] == 'HIVE Lab':
        return "Raw sequencing reads were collected as part of a fusogenicity and pathogenicity study, outlined in the following paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8942852/"
    elif row['bioproject'] ==  'PRJNA603194':
        return "SARS-CoV-2 original isolate from human lung metagenome from Wuhan."
    elif row['bioproject'] ==  'PRJEB12890':
        return "List of SRA IDs retrieved using search string in SRA 'txid211044[Organism:exp].' The SRA ID was selected based on NCBI search filters Source: RNA, Type: genome, Library layout: paired, Platform: Illumina. Project National Collection of Pathogenic Viruses (NCPV) UK sequences well-characterised, authenticated human pathogenic viruses."

dfinal=dfinal.assign(selection_notes=dfinal.apply(selection_notes, axis=1))
dfinal=dfinal.sort_values('organism_name', ascending=True)

dfinal.to_csv('ngs_id_from_api.tsv', sep = '\t', index = False)

#clean up files
os.remove('Crandallapi.json')
os.remove('hiveapi.json')
#comment out lines below to keep tsv per lab if necessary
os.remove('hiveapi.tsv')
os.remove('hiveTEMP.tsv')
os.remove('Crandallapi.tsv')
os.remove('CrandallTEMP.tsv')
