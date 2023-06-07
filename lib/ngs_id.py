import os
import requests
import json
import csv
import pandas as pd
pd.options.mode.chained_assignment = None  # removes warning message from overwriting


#Read/write tsv - writes appropriate headers for ngs_id then writes corresponding columns from ngsQC file

with open('ngsQC_HIVE.tsv', 'r') as source:
    reader = csv.reader(source, delimiter='\t')
    with open('ngs_api.tsv', 'w', newline = '') as result:
        writer=csv.writer(result, delimiter='\t')
        headings=next(reader)
        #append headers for ngs_id
        writer.writerow(['organism_name', 'intraspecific_name', 'lineage', 'genome_assembly_id', 'taxonomy_id', 'bioproject','biosample','sra_run_id','ngs_read_file_source','ref_org','selection_notes','lab_name','files_processed'])
        # These columns need to be updated/changed to reflect new columns added to ngs_id
        for r in reader:
            writer.writerow([r[0],r[1],r[2],r[3],r[4],r[10],r[11], r[12],r[14]])


df_hl = pd.read_table('ngs_api.tsv', sep='\t', on_bad_lines='skip')



df_hl=df_hl[df_hl['sra_run_id'].isnull() | ~df_hl[df_hl['sra_run_id'].notnull()].duplicated(subset='sra_run_id',keep='first')]
df_hl.lab_name = 'HIVE Lab'
df_hl.files_processed = 'ngsQC_HIVE'





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
df_hl=df_hl.assign(ref_org=df_hl.apply(ref_org, axis =1))

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
    elif row['bioproject'] == 'PRJNA550394':
        return "Selected from BioProject PRJNA550394."
    elif row['bioproject'] == 'PRJNA917703':
        return "Circulating flu strain from 2022-2023 flu season"
    elif row['bioproject'] == 'PRJNA940794':
        return "Multi-drug resistant organism with increased rate of infections beginning in late 2020."


df_hl=df_hl.assign(selection_notes=df_hl.apply(selection_notes, axis=1))
df_hl=df_hl.sort_values('organism_name', ascending=True)


#Reads table downloaded from NCBI Pathogen Isolate identifier database and removes all columns other than Biosample and isolate id
df2 = pd.read_table('ngs_id_isolates.tsv', sep='\t', on_bad_lines='skip')
df2.rename(columns={'#BioSample': 'biosample', 'Isolate': 'isolate_identifiers'}, inplace=True)

#Appends isolate id and re-orders columns
df_hl = df_hl.merge(df2, on='biosample', how='left')
df_hl = df_hl[['organism_name', 'intraspecific_name', 'lineage', 'genome_assembly_id', 'taxonomy_id', 'bioproject','biosample','sra_run_id','ngs_read_file_source','ref_org','isolate_identifiers','selection_notes','lab_name','files_processed']]



df3 = df_hl['bioproject'].drop_duplicates()



df_hl.to_csv('ngs_id_list.tsv', sep = '\t', index = False)
#List of bioproject ids - use to populate new table for isolate identifiers - then re-run
df3.to_csv('bioproject_list.tsv', index = False, header=0)

#clean up files
os.remove('ngs_api.tsv')
os.remove('bioproject_list.tsv')
