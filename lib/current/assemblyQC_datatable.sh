#!/bin/bash

#input: bash assemblyQC_datatable.sh /Users/christiewoodside/Desktop/ARGOS/mar19/all/ /Users/christiewoodside/Desktop/ARGOS/mar19/all/test_newcode.tsv
#UPDATED ON JULY 15: This update reflects the recalculations of perecent reads aligned and percent reads unaligned which was incorrect in the json outputs. This has been fixed, but some will need ot be recalcualted.
#UPDATED ON MAY 12,2025: The updates have been made to reflect how assembly is not longer pasted into the JSON. The edits allow for the genome assembly ID to be grabbed from the
#nucleotide acession now. The GCF id you see is the most updated ID accession associated with that nucleotide so it should be correct. the output table
#should be all refseq data unless there isn't any

# Created March 20, 2025
#This code is used to create the assemblyQC_* tables for ARGOSDB currently. This was because there were server issues with NCBI
read -sp "Enter your API key: " YOURAPIKEY
sleeptime_wtoken=0.14    #this will be the sleeptime set before all the API calls



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

#new as of May 9
# Clean JSON files before processing
for json_file in "$input_dir"*-qcAll.json; do
    # Remove control characters (non-printable) from JSON files and save cleaned version
    echo "Cleaning $json_file ..."
    tr -d '\000-\031' < "$json_file" > "$json_file.cleaned"
    mv "$json_file.cleaned" "$json_file"
    echo "Cleaned $json_file"
done
#-----------

# Write headers to the TSV file manually (this is specifically for assemblyQC table)
echo -e "organism_name\tinfraspecific_name\tassembled_genome_acc\tgenome_assembly_id\trepresentative_genome_acc\trepresentative_genome_org\trepresentative_genome_uniprot_acc\tlineage\ttaxonomy_id\tbco_id\tschema_version\tanalysis_platform\tanalysis_platform_object_id\tassembly_file_source\tgenomic_section\tnum_chromosomes\tnum_genes\tassembly_gc_content\tlength\tsize_gaps\tsize_contigs\tcontig_percentile\tcontig_momentum\tcoverage_contigs\tcnt_contigs\tcoverage_gaps\tcnt_gaps\tgap_percentile\tgenome_coverage\tn50\tn75\tn90\tn95\tl50\tl75\tl90\tl95\tphred_average\tcount_major_mutations\tcount_major_indels\tmutation_momentum\tindels_momentum\tmajor_mutation_momentum\tmajor_indels_momentum\talignment_anisotropy\toverhang_momentum\taligned_momentum\tentropic_momentum\treads_unaligned\tpercent_reads_unaligned\tpercent_reads_aligned\treads_aligned\tassembly_level\trpkm" > "$output_file"

# Loop through all JSON files in the provided folder (nested directory structure) matching *-qcNGS.json
for json_file in "$input_dir"*-qcAll.json; do
    echo "Grabbing values from $json_file ..."
    echo ""

    #Will be the same genome assembly ID for all rows so putting it here
    GAID=$(jq -r '.assembly // "NA"' "$json_file")
    #echo "$GAID"

    # Get total reads_aligned for this genome_assembly_id (GAID) July 25,2025
    total_aligned=$(jq -r --arg GAID "$GAID" '
    .refseq[]
    | select(.assembled_genome_acc != null)
    | .reads_aligned // 0
    ' "$json_file" | paste -sd+ - | bc)

    jq -c '.refseq[]' "$json_file" | while read -r entry; do
        # Extract assembled_genome_acc which is the nucleotide ID
        nucleotide=$(echo "$entry" | jq -r '.assembled_genome_acc // "NA"')
        #echo "$nucleotide"

        if [ "$GAID" == "NA" ]; then
            #echo "There is no genome assembly ID"
            GAID="$nucleotide"   #just making this GAID so its not as hard to go and change this everywhereeeeee
            #echo "$GAID"
        fi

        # Extract values specific to this assembled_genome_acc
        analysis_platform_object_id=$(echo "$entry" | jq -r '.analysis_platform_object_id // "NA"')
        analysis_platform=$(echo "$entry" | jq -r '.analysis_platform // "NA"')
        bco_id=$(echo "$entry" | jq -r '.id // "NA"')
        length=$(echo "$entry" | jq -r '.length // "NA"')
        size_gaps=$(echo "$entry" | jq -r '.size_gaps // "NA"')
        size_contigs=$(echo "$entry" | jq -r '.size_contigs // "NA"')
        genome_coverage=$(echo "$entry" | jq -r '.genome_coverage // "NA"')
        coverage_contigs=$(echo "$entry" | jq -r '.coverage_contigs // "NA"')
        coverage_gaps=$(echo "$entry" | jq -r '.coverage_gaps // "NA"')
        cnt_contigs=$(echo "$entry" | jq -r '.cnt_contigs // "NA"')
        cnt_gaps=$(echo "$entry" | jq -r '.cnt_gaps // "NA"')
        contig_percentile=$(echo "$entry" | jq -r '.contig_percentile // "NA"')
        gap_percentile=$(echo "$entry" | jq -r '.gap_percentile // "NA"')
        contig_momentum=$(echo "$entry" | jq -r '.contig_momentum // "NA"')

        n50=$(echo "$entry" | jq -r '.n50 // "NA"')
        l50=$(echo "$entry" | jq -r '.l50 // "NA"')
        n75=$(echo "$entry" | jq -r '.n75 // "NA"')
        l75=$(echo "$entry" | jq -r '.l75 // "NA"')
        n90=$(echo "$entry" | jq -r '.n90 // "NA"')
        l90=$(echo "$entry" | jq -r '.l90 // "NA"')
        n95=$(echo "$entry" | jq -r '.n95 // "NA"')
        l95=$(echo "$entry" | jq -r '.l95 // "NA"')

        assembly_gc_content=$(echo "$entry" | jq -r '.assembly_gc_content // "NA"')
        phred_average=$(echo "$entry" | jq -r '.phred_average // "NA"')
        count_major_mutations=$(echo "$entry" | jq -r '.count_major_mutations // "NA"')
        count_major_indels=$(echo "$entry" | jq -r '.count_major_indels // "NA"')

        mutation_momentum=$(echo "$entry" | jq -r '.mutation_momentum // "NA"')
        indels_momentum=$(echo "$entry" | jq -r '.indels_momentum // "NA"')
        major_mutation_momentum=$(echo "$entry" | jq -r '.major_mutation_momentum // "NA"')
        major_indels_momentum=$(echo "$entry" | jq -r '.major_indels_momentum // "NA"')

        alignment_anisotropy=$(echo "$entry" | jq -r '.alignment_anisotropy // "NA"')
        overhang_momentum=$(echo "$entry" | jq -r '.overhang_momentum // "NA"')
        aligned_momentum=$(echo "$entry" | jq -r '.aligned_momentum // "NA"')
        entropic_momentum=$(echo "$entry" | jq -r '.entropic_momentum // "NA"')

        reads_unaligned=$(echo "$entry" | jq -r '.reads_unaligned // "NA"')
        reads_aligned=$(echo "$entry" | jq -r '.reads_aligned // "NA"')

        #add this below July 25 to fix the values
        total_reads=$(echo "$reads_unaligned + $total_aligned" | bc)
        if [ "$total_reads" -gt 0 ]; then
            percent_reads_unaligned=$(awk -v unaligned="$reads_unaligned" -v total="$total_reads" 'BEGIN { printf "%.2f%%", (unaligned / total) * 100 }')
            percent_reads_aligned=$(awk -v aligned="$reads_aligned" -v total="$total_reads" 'BEGIN { printf "%.2f%%", (aligned / total) * 100 }')
        else
            percent_reads_unaligned="NA"
            percent_reads_aligned="NA"
        fi

        #old assignment
        # percent_reads_aligned=$(echo "$entry" | jq -r '.percent_reads_aligned // "NA"')
        # percent_reads_unaligned=$(echo "$entry" | jq -r '.percent_reads_unaligned // "NA"')
        rpkm=$(echo "$entry" | jq -r '.rpkm // "NA"')


        # Grabbing Metadata from the NCBI APIS below---------------------------------------------------------------------------------------------


        # Query the assembly database using eutils API and $GAID which is the assemblyID
        #sleep "$sleeptime_wtoken"
        SEARCH_RESULT=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=$GAID&retmode=json&api_key=$YOURAPIKEY")  #### ADD API KEY
        sleep "$sleeptime_wtoken"
        ASSEM_ID=$(echo "$SEARCH_RESULT" | jq -r '.esearchresult.idlist[0] // empty')
        sleep "$sleeptime_wtoken"
        #echo "ASSEMBLY API Response: $SEARCH_RESULT"

        ########If there is no assemblyUID associated with this nucleotide accession, lets search through nucleotide to get it
        if [ -z "$ASSEM_ID" ]; then
        #This case is usually for weird contigs like NZ_JABVAQ010000008
            echo "  No assembly UID found. Searching with nucleotide ID to find UID..."

            # If ASSEM_ID is empty, search with nucleotide
            sleep "$sleeptime_wtoken"
            SEARCH_RESULT=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=$GAID&retmode=json&api_key=$YOURAPIKEY")
            sleep "$sleeptime_wtoken"
            ASSEM_ID=$(echo "$SEARCH_RESULT" | jq -r '.esearchresult.idlist[0] // empty') #getting the nucleotide UID to add into the assembly search
            sleep "$sleeptime_wtoken"
            #echo ""
            #echo "ASSEMBLY with Nuc API Response: $SEARCH_RESULT" #getting the uid we need

            sleep "$sleeptime_wtoken"
            NUC_METADATA=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$ASSEM_ID&retmode=xml&api_key=$YOURAPIKEY")
            sleep "$sleeptime_wtoken"
            DECODED_NUC_METADATA=$(echo "$NUC_METADATA" | \
            sed 's/&lt;/</g; s/&gt;/>/g; s/&amp;/&/g' | \
            sed -E 's/(<GBQualifier_name>[^<]+<\/GBQualifier_name>)[[:space:]]*([^<]+)/\1<GBQualifier_value>\2<\/GBQualifier_value>/g' | \
            sed -E 's/<GBQualifier_value>([^<]*[<>][^<]*)<\/GBQualifier_value>/<GBQualifier_value><![CDATA[\1]]]]><![CDATA[>]]><\/GBQualifier_value>/g' | \
            sed -E 's/-->/--&gt;/g' | \
            sed -E 's/<GBFeature_location>[[:space:]]*<([0-9]+[.][.][0-9]+)<\/GBFeature_location>/<GBFeature_location>\1<\/GBFeature_location>/g' | \
            sed -E 's/<GBFeature_location>[[:space:]]*complement\(<([0-9]+[.][.][0-9]+)\)<\/GBFeature_location>/<GBFeature_location>complement(\1)<\/GBFeature_location>/g' | \
            sed -E 's/<GBFeature_location><1..\>[^<]*<\/GBFeature_location>//g' | \
            sed -E 's/<GBFeature_location>join\(<1..\>[^<]*<\/GBFeature_location>//g')

            # DECODED_NUC_METADATA=$(echo "$NUC_METADATA" | \
            # sed 's/&lt;/</g; s/&gt;/>/g; s/&amp;/&/g' | \
            # sed -E 's/(<GBQualifier_name>[^<]+<\/GBQualifier_name>)[[:space:]]*([^<]+)/\1<GBQualifier_value>\2<\/GBQualifier_value>/g')

            # Extract the entire GBXref section for debugging
            #echo "$DECODED_NUC_METADATA" | xmllint --xpath "//GBSeq/GBSeq_xrefs/GBXref" - 2>/dev/null

            # Extract Assembly ID (GCF_005280755.1) from the XML
            ASSEMBLY_ID=$(echo "$DECODED_NUC_METADATA" | xmllint --xpath "string(//GBSeq/GBSeq_xrefs/GBXref[GBXref_dbname='Assembly']/GBXref_id)" - 2>/dev/null)
            #echo "ASSEMBLY ID: $ASSEMBLY_ID"

            #finally grabbed the assembly ID from the nucleotide id. Now researching the assembly database with it
            sleep "$sleeptime_wtoken"
            SEARCH_RESULT=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=$ASSEMBLY_ID&retmode=json&api_key=$YOURAPIKEY")
            sleep "$sleeptime_wtoken"
            ASSEM_ID=$(echo "$SEARCH_RESULT" | jq -r '.esearchresult.idlist[0] // empty')
            sleep "$sleeptime_wtoken"
        
        fi

    

        #if we now have access to the correct assembly ID then
        if [[ -n "$ASSEM_ID" ]]; then

            # Query the ASSEMBLY metadata to get more information to fill out the table
            sleep "$sleeptime_wtoken"
            ASSEM_METADATA=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=$ASSEM_ID&retmode=xml&api_key=$YOURAPIKEY")
            sleep "$sleeptime_wtoken"
            
            # Decode HTML entities in the ASSEM_METADATA
            DECODED_ASSEM_METADATA=$(echo "$ASSEM_METADATA" | sed 's/&lt;/</g; s/&gt;/>/g; s/&amp;/&/g')
            #echo "$DECODED_ASSEM_METADATA" >> assembly_metadata.txt

            # Extract LastMajorReleaseAccession (the GCF ID) to be in the table as assembly ID
            GCF_ID=$(echo "$DECODED_ASSEM_METADATA" \
            | xmllint --xpath 'string(//*[local-name()="LastMajorReleaseAccession"])' - 2>/dev/null)
            #echo "$GCF_ID"
            TAXID=$(echo "$DECODED_ASSEM_METADATA" | xmllint --xpath "//*[local-name()='Taxid']/text()" - 2>/dev/null)
            ORGANISM=$(echo "$DECODED_ASSEM_METADATA" | xmllint --xpath "//*[local-name()='SpeciesName']/text()" - 2>/dev/null)
            ASSEMBLY_STATUS=$(echo "$DECODED_ASSEM_METADATA" | xmllint --xpath "//*[local-name()='AssemblyStatus']/text()" - 2>/dev/null)
            

            # Extract the content of the Meta tag (which includes CDATA)
            META_CONTENT=$(echo "$DECODED_ASSEM_METADATA" | xmllint --xpath "//*[local-name()='Meta']/text()" - 2>/dev/null)
            # Remove CDATA wrapper by replacing the CDATA tags with their content
            CLEAN_META_CONTENT=$(echo "$META_CONTENT" | sed 's/<!\[CDATA\[\(.*\)\]\]>/\1/')
            # Now parse the cleaned META content
            # Ensure the content is valid XML by wrapping it in a proper root element
            CLEAN_META_XML="<root>$CLEAN_META_CONTENT</root>"
            
            # Extract chromosome count from the cleaned and wrapped META content
            CHROMOSOME_COUNT=$(echo "$CLEAN_META_XML" | xmllint --xpath 'string(//Stat[@category="chromosome_count" and @sequence_tag="all"])' -)
            #echo "chromosome counts: $CHROMOSOME_COUNT"
            
        fi




        # Query the nucleotide database using eutils API and $nucleotide which is the nucleotide assecion iD
        sleep "$sleeptime_wtoken"
        SEARCH_RESULT=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=$nucleotide&retmode=json&api_key=$YOURAPIKEY")
        sleep "$sleeptime_wtoken"
        NUC_ID=$(echo "$SEARCH_RESULT" | jq -r '.esearchresult.idlist[0] // empty')
        sleep "$sleeptime_wtoken"
        #echo "NUC API Response: $SEARCH_RESULT"

        if [[ -n "$NUC_ID" ]]; then
            # Query the NUCLEOTIDE metadata to get more information to fill out the table
            sleep "$sleeptime_wtoken"
            NUC_METADATA=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$NUC_ID&retmode=xml&api_key=$YOURAPIKEY")
            sleep "$sleeptime_wtoken"
            #echo "$NUC_METADATA"

            # Decode HTML entities in the NUC_METADATA -> some tags/f;ags were causing parsing errors which is why these are so long
            #        - Lines commented out worked, but there were some instances where the aprsing wouldn't work
            #DECODED_NUC_METADATA=$(echo "$NUC_METADATA" | sed 's/&lt;/</g; s/&gt;/>/g; s/&amp;/&/g' | sed 's/<GBFeature_location><[0-9]*\.\.[0-9]*<\/GBFeature_location>//g')
                    # - Second best working one:
            #DECODED_NUC_METADATA=$(echo "$NUC_METADATA" | sed 's/&lt;/</g; s/&gt;/>/g; s/&amp;/&/g' | sed -e 's/<GBFeature_location><[0-9]*\.\.[0-9]*<\/GBFeature_location>//g' -e 's/<GBFeature_location>complement(<[^<]*<\/GBFeature_location>//g' -e 's/<GBQualifier_value>.*< J Loper.*<\/GBQualifier_value>//g')
                    # - Current best and working:
            # DECODED_NUC_METADATA=$(echo "$NUC_METADATA" | \
            # sed 's/&lt;/</g; s/&gt;/>/g; s/&amp;/&/g' | \
            # sed -E 's/(<GBQualifier_name>[^<]+<\/GBQualifier_name>)[[:space:]]*([^<]+)/\1<GBQualifier_value>\2<\/GBQualifier_value>/g' | \
            # sed -E 's/<GBQualifier_value>([^<]*[<>][^<]*)<\/GBQualifier_value>/<GBQualifier_value><![CDATA[\1]]]]><![CDATA[>]]><\/GBQualifier_value>/g' | \
            # sed -E 's/-->/--&gt;/g' | \
            # sed -E 's/<GBFeature_location>[[:space:]]*<([0-9]+[.][.][0-9]+)<\/GBFeature_location>/<GBFeature_location>\1<\/GBFeature_location>/g' | \
            # sed -E 's/<GBFeature_location>[[:space:]]*complement\(<([0-9]+[.][.][0-9]+)\)<\/GBFeature_location>/<GBFeature_location>complement(\1)<\/GBFeature_location>/g' | \
            # sed -E 's/<GBFeature_location><([0-9]+[.][.][0-9]+)<\/GBFeature_location>/<GBFeature_location>\1<\/GBFeature_location>/g' | \
            # sed -E 's/<GBFeature_location>[[:space:]]*<([0-9]+[.][.][0-9]+)<\/GBFeature_location>/<GBFeature_location>\1<\/GBFeature_location>/g')

            DECODED_NUC_METADATA=$(echo "$NUC_METADATA" | \
            sed 's/&lt;/</g; s/&gt;/>/g; s/&amp;/&/g' | \
            sed -E 's/(<GBQualifier_name>[^<]+<\/GBQualifier_name>)[[:space:]]*([^<]+)/\1<GBQualifier_value>\2<\/GBQualifier_value>/g' | \
            sed -E 's/<GBQualifier_value>([^<]*[<>][^<]*)<\/GBQualifier_value>/<GBQualifier_value><![CDATA[\1]]]]><![CDATA[>]]><\/GBQualifier_value>/g' | \
            sed -E 's/-->/--&gt;/g' | \
            sed -E 's/<GBFeature_location>[[:space:]]*<([0-9]+[.][.][0-9]+)<\/GBFeature_location>/<GBFeature_location>\1<\/GBFeature_location>/g' | \
            sed -E 's/<GBFeature_location>[[:space:]]*complement\(<([0-9]+[.][.][0-9]+)\)<\/GBFeature_location>/<GBFeature_location>complement(\1)<\/GBFeature_location>/g' | \
            sed -E 's/<GBFeature_location><1..\>[^<]*<\/GBFeature_location>//g' | \
            sed -E 's/<GBFeature_location>join\(<1..\>[^<]*<\/GBFeature_location>//g')


            #gives error outputs
            #echo "$DECODED_NUC_METADATA" | xmllint --noout -
            #echo "$DECODED_NUC_METADATA" > debug_nucleotide.xml
            xmllint --format debug_nucleotide.xml > formatted_output.xml
            #echo "$DECODED_NUC_METADATA" >> nucleotide_metadata.txt


            # Extract values to add to the table--------
            #LINEAGE=$(echo "$DECODED_NUC_METADATA" | xmllint --xpath "string(//GBSeq_taxonomy)" - 2>/dev/null)
                #This one is more direct/specific than the above one
            LINEAGE=$(echo "$DECODED_NUC_METADATA" | xmllint --xpath "string(//GBSeq/GBSeq_taxonomy)" - 2>/dev/null)
            GENES_TOTAL=$(echo "$DECODED_NUC_METADATA" | grep -o '; Genes (total) :: [0-9,]\+' | sed 's/; Genes (total) :: //')
            GENE_SEC=$(echo "$DECODED_NUC_METADATA" | xmllint --xpath 'string(//GBSeq_definition)' - 2>/dev/null)
        



            #How I am getting the genomic section value for the table

            OUTPUT=""
            if [[ -z "$GENE_SEC" ]]; then
                OUTPUT="  "
                echo "***** No GBSEq_Definition (for Genome_Section) found *****"

            elif [[ "$GENE_SEC" == *"Influenza A virus"* && "$GENE_SEC" == *"))"* ]]; then
                # Extract the infraspecific name — what's inside the outermost parentheses after "Influenza A virus"
                INFRA_NAME=$(echo "$GENE_SEC" | sed -n 's/.*Influenza A virus (\([^)]*\)).*/\1/p')
                INFRA_NAME="${INFRA_NAME}))"

                # Use awk to get the part after the last )) and before the last comma
                TRIMMED=$(echo "$GENE_SEC" | awk -F'\\)\\)' '{print $NF}' | awk -F',' 'NF{NF--; print}' OFS=,)
                #TRIMMED=$(echo "$GENE_SEC" | awk -F'\\)' '{print $NF}' | awk -F',' 'NF{NF--; print}' OFS=,)
                OUTPUT=$(echo "$TRIMMED" | xargs)  # Trim leading/trailing whitespace
                echo "-------- Influenza Segment Extracted: $OUTPUT"


            elif [[ $GENE_SEC != *","* ]]; then
                echo "-------- Genomic_Section w/o ,: $GENE_SEC"

                if [[ $GENE_SEC == *"chromosome"* ]]; then
                    OUTPUT="chromosome"
                
                elif [[ $GENE_SEC == *"plasmid unnamed"* ]]; then
                    # Extract the number from "plasmid unnamedX"
                    if [[ $GENE_SEC =~ plasmid\ unnamed([0-9]+) ]]; then
                        OUTPUT="plasmid unnamed${BASH_REMATCH[1]}"
                    else
                        OUTPUT="plasmid unnamed"
                    fi
                
                elif [[ $GENE_SEC == *"plasmid"* && $GENE_SEC != *"plasmid unnamed"* ]]; then
                    OUTPUT="plasmid"
                
                elif [[ $GENE_SEC == *"unnamed"* ]] && [[ $GENE_SEC != *"plasmid unnamed"* ]]; then
                    # Extract the number from "unnamedX"
                    if [[ $GENE_SEC =~ unnamed([0-9]+) ]]; then
                        OUTPUT="unnamed${BASH_REMATCH[1]}"
                    else
                        OUTPUT="unnamed"
                    fi
                
                else
                    OUTPUT=""
                fi

            #if there is a , in the GBSeq_definition name (which there usually is for FDA-ARGOS strains)
            else
                RESULT=$(echo $GENE_SEC | cut -d ',' -f1 | xargs)  # Extract first part before comma
                #echo "Genomic_Section before , result: $RESULT"
                echo "        $RESULT"

                if [[ $GENE_SEC == *"contig"* || $GENE_SEC == *"unitig_0_quiver_pilon"* || \
                    $GENE_SEC == *"trim_quiver_pilon"* || $GENE_SEC == *"_quiver_pilon"* || \
                    $GENE_SEC == *"_arrow_pilon"* ]]; then
                    OUTPUT="contig"

                elif [[ $GENE_SEC == *"plasmid unnamed"* ]]; then
                    # Extract the number from "plasmid unnamedX"
                    if [[ $GENE_SEC =~ plasmid\ unnamed([0-9]+) ]]; then
                        OUTPUT="plasmid unnamed${BASH_REMATCH[1]}"
                    else
                        OUTPUT="plasmid unnamed"
                    fi

                elif [[ $GENE_SEC == *"chromosome"* ]]; then
                    OUTPUT="chromosome"

                elif [[ $GENE_SEC == *"plasmid"* && $GENE_SEC != *"plasmid unnamed"* ]]; then
                    OUTPUT="plasmid"

                elif [[ $GENE_SEC == *"unnamed"* ]] && [[ $GENE_SEC != *"plasmid unnamed"* ]]; then
                    # Extract the number from "unnamedX"
                    if [[ $GENE_SEC =~ unnamed([0-9]+) ]]; then
                        OUTPUT="unnamed${BASH_REMATCH[1]}"
                    else
                        OUTPUT="unnamed"
                    fi

                else
                    OUTPUT="$RESULT"
                fi
            fi
            # echo "Genomic Section Output: $OUTPUT"
            #echo ""

        fi

         # Append extracted values to the output TSV file
        echo -e "$ORGANISM\t$INFRA_NAME\t$nucleotide\t$GCF_ID\t\t\t\t$LINEAGE\t$TAXID\t\tv1.6\t$analysis_platform\t$analysis_platform_object_id\tNCBI\t$OUTPUT\t$CHROMOSOME_COUNT\t$GENES_TOTAL\t$assembly_gc_content\t$length\t$size_gaps\t$size_contigs\t$contig_percentile\t$contig_momentum\t$coverage_contigs\t$cnt_contigs\t$coverage_gaps\t$cnt_gaps\t$gap_percentile\t$genome_coverage\t$n50\t$n75\t$n90\t$n95\t$l50\t$l75\t$l90\t$l95\t$phred_average\t$count_major_mutations\t$count_major_indels\t$mutation_momentum\t$indels_momentum\t$major_mutation_momentum\t$major_indels_momentum\t$alignment_anisotropy\t$overhang_momentum\t$aligned_momentum\t$entropic_momentum\t$reads_unaligned\t$percent_reads_unaligned\t$percent_reads_aligned\t$reads_aligned\t$ASSEMBLY_STATUS\t$rpkm" >> "$output_file"
    done
done

# Indicate script completion
echo ""
echo "Processing complete. Output written to $output_file"
