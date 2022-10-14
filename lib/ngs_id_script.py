import csv
import pandas as pd
pd.options.mode.chained_assignment = None  # removes warning message from overwriting --- I think

#Read through HIVE Lab ngs_QC and pull out relevant data
with open("c:/users/jgerg/data.argosdb/home/template_ngsQC_HL.tsv", 'r') as inFile1:
    tsvreader = csv.reader(inFile1, delimiter="\t")
    with open('trimmed_HL.tsv', 'w', newline='') as outFile:
        tsvwriter = csv.writer(outFile, delimiter="\t")
        #remove headers from ngs_QC file
        headings=next(tsvreader)
        #append headers for ngs_id
        tsvwriter.writerow(['organism_name', 'leaf_node', 'genome_assembly_id', 'taxonomy_id', 'bioproject','sra_run_id','ref_org','selection_notes','lab_name','files_processed'])
        for row in tsvreader:
            tsvwriter.writerow([row[0],row[1],row[2],row[3],row[9],row[11]])
##
##
### Removes duplicate genome assembly id's but keeps any row that does not have an assembly id
### Need to add  this code - it doesnt correctly remove duplicates of samples with no genome assembly id
### For example the different SARS variants - none of them have genome assembly ids - maybe also removing duplicate SRRs will suffice
###
###
data_hl = pd.read_table('trimmed_HL.tsv', sep='\t')
df_hl = data_hl
df_hl=df_hl[df_hl['genome_assembly_id'].isnull() | ~df_hl[df_hl['genome_assembly_id'].notnull()].duplicated(subset='genome_assembly_id',keep='first')]
df_hl.lab_name = 'Hive Lab'
df_hl.files_processed = 'ngsQC_HL'

###Will print tsv for unique genome assemly ids for HIVE Lab
#df_hl.to_csv('trimmed_HL.tsv', sep = '\t', index = False)

#Repeat for Pond Lab
with open("c:/users/jgerg/data.argosdb/home/template_ngsQC_Pond.tsv", 'r') as inFile2:
    tsvreader = csv.reader(inFile2, delimiter="\t")
    with open('trimmed_Pond.tsv', 'w', newline='') as outFile:
        tsvwriter = csv.writer(outFile, delimiter="\t")
        #remove headers from ngs_QC file
        headings=next(tsvreader)
        #append headers for ngs_id
        tsvwriter.writerow(['organism_name', 'leaf_node', 'genome_assembly_id', 'taxonomy_id', 'bioproject','sra_run_id','ref_org','selection_notes','lab_name','files_processed'])
        for row in tsvreader:
            tsvwriter.writerow([row[0],row[1],row[2],row[3],row[9],row[11]])
data_p = pd.read_table('trimmed_Pond.tsv', sep='\t')
df_p = data_p
df_p=df_p[df_p['genome_assembly_id'].isnull() | ~df_p[df_p['genome_assembly_id'].notnull()].duplicated(subset='genome_assembly_id',keep='first')]
df_p.lab_name = 'Pond Lab'
df_p.files_processed = 'ngsQC_Pond'

###Will print tsv for unique genome assemly ids for Pond Lab
#df_p.to_csv('trimmed_Pond.tsv', sep = '\t', index = False)


#Repeat for Crandall Lab
with open("c:/users/jgerg/data.argosdb/home/template_ngsQC_Crandall.tsv", 'r') as inFile3:
    tsvreader = csv.reader(inFile3, delimiter="\t")
    with open('trimmed_Crandall.tsv', 'w', newline='') as outFile:
        tsvwriter = csv.writer(outFile, delimiter="\t")
        #remove headers from ngs_QC file
        headings=next(tsvreader)
        #append headers for ngs_id
        tsvwriter.writerow(['organism_name', 'leaf_node', 'genome_assembly_id', 'taxonomy_id', 'bioproject','sra_run_id','ref_org','selection_notes','lab_name','files_processed'])
        for row in tsvreader:
            tsvwriter.writerow([row[0],row[1],row[2],row[3],row[9],row[11]])
data_c = pd.read_table('trimmed_Crandall.tsv', sep='\t')
df_c = data_c
df_c=df_c[df_c['genome_assembly_id'].isnull() | ~df_c[df_c['genome_assembly_id'].notnull()].duplicated(subset='genome_assembly_id',keep='first')]
df_c.lab_name = 'Crandall Lab'
df_c.files_processed = 'ngsQC_Crandall'

###Will print tsv for unique genome assemly ids for HIVE Lab
#df_c.to_csv('trimmed_Crandall.tsv', sep = '\t', index = False)



dfs = [df_hl, df_p, df_c]
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
    if row['bioproject'] == 'PRJNA231221': #ARGOS Bioproject
           return "Belongs to FDA-ARGOS PRJNA231221."
    ## add ifelse to populate selection_notes for SARS

#Drop bioproject column - was only used to populate selection_notes - isnt part of ngs_id_list
dfinal=dfinal.assign(selection_notes=dfinal.apply(selection_notes, axis=1))
dfinal=dfinal.drop('bioproject', axis=1)
dfinal=dfinal[dfinal['genome_assembly_id'].isnull() | ~dfinal[dfinal['genome_assembly_id'].notnull()].duplicated(subset='genome_assembly_id',keep='first')]
dfinal.to_csv('c:/users/jgerg/data.argosdb/home/test_ngs_id_list.tsv', sep = '\t', index = False)
