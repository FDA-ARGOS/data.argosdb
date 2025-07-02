#!/bin/bash

# Pull all read information associated with a list of assemblies (SRR, ERR, and DRR).
# Note that BioSample is commented out, because this will be pulled by HIVE automatically.
# This script is written for anonymous access (sleep .8 seconds per query) since it's a public repo (and I haven't set up .secrets yet). Use your NCBI API auth key, if you like.

while read -u3 line
do

    echo "Reading $line..."

	sra="$(esearch -db bioproject -query $line | elink -target sra | esummary)" || echo "$line failed."
    sleep .4

	srr="$(echo "${sra}" | grep -E -o 'SRR[0-9]{4,}')"
    if [ -z "${srr}" ];then
        echo "Nothing found for $line."
    else
        echo "${srr}" >> readList2.txt
    fi

	err="$(echo "${sra}" | grep -E -o 'SRR[0-9]{4,}')"
    if [ -z "${err}" ];then
        echo "Nothing found for $line."
    else
        echo "${err}" >> readList2.txt
    fi

	drr="$(echo "${sra}" | grep -E -o 'SRR[0-9]{4,}')"
    if [ -z "${drr}" ];then
        echo "Nothing found for $line."
    else
        echo "${drr}" >> readList2.txt
    fi

	sleep 1

done 3< rajaBioprojects.txt

while read -u3 row
do

    echo "Reading $row..."

	assembly="$(esearch -db sra -query $row | elink -target assembly | esummary)" | grep "FtpPath_RefSeq" | sed -r 's/.*>(ftp:\/\/ftp.ncbi.nlm.nih.gov\/genomes\/all\/.+.+)<.*/\1/' | rev | cut -d'/' -f 1 | rev >> assemblies.txt || echo "$row failed."
    sleep .4

done 3< readList2.txt