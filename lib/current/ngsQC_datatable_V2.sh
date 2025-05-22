#!/bin/bash

# input example: bash ngsQC_datatable_V2.sh /Users/USER/Desktop/ARGOS/feb5/ngs/ /Users/USER/Desktop/ARGOS/feb5/ngs/newcode.tsv
# Updated may 20, 2025 by Christie Woodside
#This code is used to create the ngsQC_HIVE3 table for ARGOSDB currently. This was because there were server issues with NCBI


#bash ngsQC_datatable_V2.sh /Users/christiewoodside/Desktop/ARGOS/may21/ngs/ /Users/christiewoodside/Desktop/ARGOS/may21/ngs/ngs_may21_test2.tsv

# Check if jq is installed
if ! command -v jq &> /dev/null; then
    echo "jq is required but not installed. Install jq and try again."
    exit 1
fi

# Input directory and output file (output file is optional, default will be used if not provided)
input_dir="$1"
output_file="${2:-'combined_data.tsv'}"  # Use provided output file, or default to 'ncbi_data.tsv'

# Check if the input directory is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <input_directory> [<output_file>]"
    exit 1
fi

# Check if the input directory exists
if [ ! -d "$input_dir" ]; then
    echo "Input directory $input_dir does not exist. Exiting."
    exit 1
fi

# Ensure output file is empty before writing
> "$output_file"

#---------------------------------------------------------------------------------------------------------------




#_-------_-------_---------_---------_---------_--------_-------_--------_--------_---------_---------_--------_---

# Write headers to the TSV file manually (can update later if needed)
echo -e "organism_name\tinfraspecific_name\tlineage\trepresentative_genome_org\trepresentative_genome_acc\trepresentative_genome_uniprot_acc\tgenome_assembly_id\ttaxonomy_id\tbco_id\tschema_version\tanalysis_platform\tanalysis_platform_object_id\tstrain\tbioproject\tbiosample\tsra_run_id\tngs_read_file_name\tngs_read_file_source\tngs_gc_content\tavg_phred_score\tmin_read_length\tnum_reads\tnum_reads_unique\tavg_read_length\tmax_read_length\tmax_duplicate_read\tstrategy\tinstrument\tid_method\tcoding_system\tpercent_coding\tpercent_not_coding\tcomplexity_percent\tnon_complex_percent\tstdev_quality\tavg_quality_a\tavg_quality_t\tavg_quality_g\tavg_quality_c\tavg_quality_n\tcount_a\tcount_c\tcount_g\tcount_t\tcount_n\tpercent_a\tpercent_c\tpercent_g\tpercent_t\tpercent_n\tcount_all_WN\tcount_all" > "$output_file"

# Loop through all JSON files in the provided folder (nested directory structure) matching *-qcNGS.json
for json_file in "$input_dir"*-qcNGS.json; do
    echo "Grabbing values from $json_file ..."

    # Refactor: collapse per-field jq calls into one, keep assembled_genome_acc and SRR_ID derivation
    jq -r '
        .ngsqc[]
        | [
            .assembled_genome_acc // "NA",
            .analysis_platform_object_id // "NA",
            .ngs_gc_content // "NA",
            .avg_phred_score // "NA",
            .min_read_length // "NA",
            .num_reads // "NA",
            .num_reads_unique // "NA",
            .avg_read_length // "NA",
            .max_read_length // "NA",
            .max_duplicate_read // "NA",
            .coding_system // "NA",
            .percent_coding // "NA",
            .percent_non_coding // "NA",
            .complexity_percent // "NA",
            .non_complexity_percent // "NA",
            .stdev_quality // "NA",
            .bases.avg_quality_a // "NA",
            .bases.avg_quality_t // "NA",
            .bases.avg_quality_g // "NA",
            .bases.avg_quality_c // "NA",
            .bases.avg_quality_n // "NA",
            .bases.count_a // "NA",
            .bases.count_c // "NA",
            .bases.count_g // "NA",
            .bases.count_t // "NA",
            .bases.count_n // "NA",
            .bases.percent_a // "NA",
            .bases.percent_c // "NA",
            .bases.percent_g // "NA",
            .bases.percent_t // "NA",
            .bases.percent_n // "NA",
            .count_all_WN // "NA",
            .count_all        // "NA"
            ]
        | @tsv
    ' "$json_file" \
    | while IFS=$'\t' read -r \
        ngs_read_file_name \
        obj_id \
        ngs_gc_content \
        avg_phred_score \
        min_read_length \
        num_reads \
        num_reads_unique \
        avg_read_length \
        max_read_length \
        max_duplicate_read \
        coding_system \
        percent_coding \
        percent_not_coding \
        complexity_percent \
        non_complex_percent \
        stdev_quality \
        avg_quality_a \
        avg_quality_t \
        avg_quality_g \
        avg_quality_c \
        avg_quality_n \
        count_a \
        count_c \
        count_g \
        count_t \
        count_n \
        percent_a \
        percent_c \
        percent_g \
        percent_t \
        percent_n \
        count_all_WN \
        count_all; do

    # Derive SRR_ID from the assembled_genome_acc
    SRR_ID="${ngs_read_file_name%%_*}"
    #echo "DEBUG: extracted SRR_ID=[$SRR_ID]"
    


        SEARCH_RESULT=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=$SRR_ID&retmode=json&api_key=bfbde99c962d228023e8d62a078bdb12d108")
        SRA_ID=$(echo "$SEARCH_RESULT" | jq -r '.esearchresult.idlist[0] // empty')
        #echo "$SRA_ID"
        #echo "SRA API Response: $SEARCH_RESULT"



        if [[ -n "$SRA_ID" ]]; then
            
            # Query the SRA metadata to get more information about the SRA ID
            SRA_METADATA=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=$SRA_ID&retmode=xml&api_key=bfbde99c962d228023e8d62a078bdb12d108")
            
            # Decode HTML entities in the SRA_METADATA
            DECODED_SRA_METADATA=$(echo "$SRA_METADATA" | sed 's/&lt;/</g; s/&gt;/>/g; s/&amp;/&/g')
            #echo "$DECODED_SRA_METADATA"

            # Extract relevant details from the metadata
            ORGANISM=$(echo "$DECODED_SRA_METADATA" | xmllint --xpath "string(//DocSum/Item[@Name='ExpXml']/descendant::Organism/@ScientificName)" - 2>/dev/null)
            ORGANISM_TAXID=$(echo "$DECODED_SRA_METADATA" | xmllint --xpath "string(//DocSum/Item[@Name='ExpXml']/descendant::Organism/@taxid)" - 2>/dev/null)
            BIOSAMPLE_ID=$(echo "$DECODED_SRA_METADATA" | xmllint --xpath "string(//DocSum/Item[@Name='ExpXml']/descendant::Biosample)" - 2>/dev/null)
            INSTRUMENT=$(echo "$DECODED_SRA_METADATA" | xmllint --xpath "string(//DocSum/Item[@Name='ExpXml']/descendant::Platform/@instrument_model)" - 2>/dev/null)
            STRATEGY=$(echo "$DECODED_SRA_METADATA" | xmllint --xpath "string(//DocSum/Item[@Name='ExpXml']/descendant::Library_descriptor/descendant::LIBRARY_STRATEGY)" - 2>/dev/null)
            SOURCE=$(echo "$DECODED_SRA_METADATA" | xmllint --xpath "string(//DocSum/Item[@Name='ExpXml']/descendant::Library_descriptor/descendant::LIBRARY_SOURCE)" - 2>/dev/null)
            BIOPROJECT=$(echo "$DECODED_SRA_METADATA" | xmllint --xpath "string(//DocSum/Item[@Name='ExpXml']/descendant::Bioproject)" - 2>/dev/null)



            # this gets me the lineage for the organisms 
            ENCODED_ORGANISM=$(echo "$ORGANISM" | sed 's/ /%20/g')

            tax_META=$(curl -s "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/$ENCODED_ORGANISM")
            LINEAGE=$(echo "$tax_META" | jq -r '.[0].lineage' | sed 's/; $//')
            #echo "$LINEAGE"


            # Query the Biosample database to get more details using the Biosample ID
            if [[ -n "$BIOSAMPLE_ID" ]]; then

                #Getting Biosample metadata that is needed
                BIOSAMPLE_SEARCH_RESULT=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=biosample&term=$BIOSAMPLE_ID&retmode=json&api_key=bfbde99c962d228023e8d62a078bdb12d108")
                
                # Extract the Biosample UID from the search result
                BIOSAMPLE_UID=$(echo "$BIOSAMPLE_SEARCH_RESULT" | jq -r '.esearchresult.idlist[0] // empty')

                if [[ -n "$BIOSAMPLE_UID" ]]; then
                    BIOSAMPLE_METADATA=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=biosample&id=$BIOSAMPLE_UID&retmode=xml&api_key=bfbde99c962d228023e8d62a078bdb12d108")
                    BIOSAMPLE_METADATA_CLEAN=$(echo "$BIOSAMPLE_METADATA" | sed 's/<!DOCTYPE[^>]*>//g')
                    #echo "$BIOSAMPLE_METADATA_CLEAN"
                    STRAIN=$(echo "$BIOSAMPLE_METADATA_CLEAN" | xmlstarlet sel -t -v "//Attribute[@attribute_name='strain']" -n)

                    if [[ -z "$STRAIN" ]]; then
                        STRAIN=$(echo "$BIOSAMPLE_METADATA_CLEAN" | xmlstarlet sel -t -v "//Infraspecies[contains(., 'strain:')]" | sed 's/strain: //')
                    fi
                
                    IDENTIFICATION_METHOD=$(python3 -c "
import time
from Bio import Entrez
import xmltodict

Entrez.email = 'christie.woodside@email.gwu.edu'

def bsMeta(bs_term, sleeptime):
    ''' gets additional biosample information to add to the tsv'''
    if not bs_term:
        print(f'Error: No valid biosample term provided.')
        return []
    search = Entrez.esearch(db='biosample', term=bs_term, retmode='xml')
    time.sleep(sleeptime)
    record = Entrez.read(search)
    time.sleep(sleeptime)
    
    if 'IdList' not in record or not record['IdList']:
        print(f'No biosample found for term: {bs_term}')
        return []

    bs_id = record['IdList'][0]
    
    info = Entrez.esummary(db='biosample', id=bs_id)
    time.sleep(sleeptime)
    record = Entrez.read(info)
    time.sleep(sleeptime)
    r = record['DocumentSummarySet']['DocumentSummary'][0]
    
    sd = r['SampleData']
    sd_json = xmltodict.parse(sd)['BioSample']
    
    idm_id = ''
    for id in sd_json['Attributes']['Attribute']:
        if id.get('@attribute_name') == 'identification method':
            idm_id = id['#text']
    
    return idm_id
    
print(bsMeta('$BIOSAMPLE_ID', 0.75))")

                        IDENTIFICATION_METHOD=$(echo "$IDENTIFICATION_METHOD" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

                else
                    STRAIN=""
                fi
            else
                STRAIN=""
            fi


            # ---- START: replace GAID pull here ----
            # 1) BioSample → Assembly
            if [[ -n "$BIOSAMPLE_UID" ]]; then
                #echo "$BIOSAMPLE_UID"
                BS_LINKS=$(curl -s \
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=biosample&db=assembly&id=$BIOSAMPLE_UID&retmode=json&api_key=bfbde99c962d228023e8d62a078bdb12d108")

                #echo "DEBUG raw BS_LINKS →"
                #echo "$BS_LINKS" | jq .
                ASSEMBLY_UID=$(echo "$BS_LINKS" \
                    | jq -r '.linksets[0]?.linksetdbs[0]?.links[0] // empty')
                #echo "$ASSEMBLY_UID"    
            fi

            # 2) fallback SRA → Assembly
            if [[ -z "$ASSEMBLY_UID" && -n "$SRA_ID" ]]; then
                SR_LINKS=$(curl -s \
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=sra&db=assembly&id=$SRA_ID&retmode=json&api_key=bfbde99c962d228023e8d62a078bdb12d108")
                ASSEMBLY_UID=$(echo "$SR_LINKS" \
                    | jq -r '.linksets[0].linksetdbs[0].links[0] // empty')
            fi

            # 3) fetch the accession
            if [[ -n "$ASSEMBLY_UID" ]]; then
                SUMMARY=$(curl -s \
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=$ASSEMBLY_UID&retmode=json&api_key=bfbde99c962d228023e8d62a078bdb12d108")
                GAID=$(echo "$SUMMARY" \
                    | jq -r --arg id "$ASSEMBLY_UID" \
                    '.result[$id].assemblyaccession // .result[$id].AssemblyAccession // "NA"')
            else
                GAID="NA"
            fi
            echo "DEBUG: final GAID → $GAID"
            # ---- END: replace GAID pull here ----




            # Save the result to the output file, adding the schema_version column
            echo -e "$ORGANISM\t\t$LINEAGE\t\t\t\t$GAID\t$ORGANISM_TAXID\t\tv1.6\tHIVE\t$obj_id\t$STRAIN\t$BIOPROJECT\t$BIOSAMPLE_ID\t$SRR_ID\t$ngs_read_file_name\tSRA\t$ngs_gc_content\t$avg_phred_score\t$min_read_length\t$num_reads\t$num_reads_unique\t$avg_read_length\t$max_read_length\t$max_duplicate_read\t$STRATEGY\t$INSTRUMENT\t$IDENTIFICATION_METHOD\t$coding_system\t$percent_coding\t$percent_not_coding\t$complexity_percent\t$non_complex_percent\t$stdev_quality\t$avg_quality_a\t$avg_quality_t\t$avg_quality_g\t$avg_quality_c\t$avg_quality_n\t$count_a\t$count_c\t$count_g\t$count_t\t$count_n\t$percent_a\t$percent_c\t$percent_g\t$percent_t\t$percent_n\t$count_all_WN\t$count_all" >> "$output_file"
        else
            echo "         No matching SRA found for ShortRead: $SRR_ID"
            echo -e "$ORGANISM\tNo matching SRA found\t$LINEAGE\t\t\t\t$GAID\t$ORGANISM_TAXID\t\tv1.6\tHIVE\t$obj_id\t$STRAIN\t$BIOPROJECT\t$BIOSAMPLE_ID\t\t$ngs_read_file_name\tSRA" >> "$output_file"
        fi
    done
    
done


echo "Data processing completed. Output saved to $output_file."
