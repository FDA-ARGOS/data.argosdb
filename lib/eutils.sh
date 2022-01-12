# SRX -> SRA_NGSQC
esearch -db sra -query SRX1382142 | efetch -format xml 





# # #Get all BioSamples that have an assembly
# esearch -db bioproject -query PRJNA231221 | \
# efetch -format xml | \
#   xtract -pattern ProjectDescr -element Title -deq "\n\t" \
#     -group LocusTagPrefix -if @assembly_id  -ASS @assembly_id -SAMP @biosample_id \
#       -element "&ASS" "&SAMP" LocusTagPrefix -deq "\n\t" | \
#         esearch -db biosample -query "&SAMP" | \
#         elink -target sra | efetch -format xml | \
#           xtract -pattern EXPERIMENT_PACKAGE -element EXPERIMENT@alias SAMPLE@accession EXPERIMENT@accession

# # Get all SRA runs for a given biosample
# esearch -db biosample -query SAMN03255469 | \
# elink -target sra |  efetch -format xml | \
#     xtract -pattern EXPERIMENT_PACKAGE -element EXPERIMENT@alias SAMPLE@accession EXPERIMENT@accession

# efetch -format xml | \
# xtract -pattern ProjectDescr  -block LocusTagPrefix\
#    -element @assembly_id @biosample_id LocusTagPrefix


# Get all SRA runs for a Biosample based on a BioSample ID
esearch -db sra -query "SAMN24370922" | \
efetch -format docsum | \
xtract -pattern Runs -ACC @acc  -element "&ACC"

# esearch -db sra -query SRR17325261 | \
# efetch -format xml | \
# xtract -pattern EXPERIMENT_PACKAGE -block SAMPLE -element -KEY PRIMARY_ID DB ID LABEL

# efetch -db nuccore -id NM_000084.3 -format gbc |
# xtract -pattern INSDSeq -element INSDSeq_accession-version \
#   -group INSDFeature -KEY INSDFeature_key \
#     -block INSDQualifier -deq "\n\t" \
#       -element "&KEY" INSDQualifier_name INSDQualifier_value