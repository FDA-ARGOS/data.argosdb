#!/bin/bash
echo "Assembly;SRR;ERR;DRR;BioSample"
while read -u3 line || [ "$line" ]
do

        assembly="$(esearch -db bioproject -query $line | elink -target assembly | esummary | grep "FtpPath_RefSeq" | sed -r 's/.*>(ftp:\/\/ftp.ncbi.nlm.nih.gov\/genomes\/all\/.+.+)<.*/\1/' | rev | cut -d'/' -f 1 | rev | tr '\n' ',' | sed 's/.$//')"

        sra="$(esearch -db bioproject -query $line | elink -target sra | esummary)"
        srr="$(echo "${sra}" | grep -E -o 'SRR[0-9]{4,}' | tr '\n' ',' | sed 's/.$//')"
        err="$(echo "${sra}" | grep -E -o 'ERR[0-9]{4,}' | tr '\n' ',' | sed 's/.$//')"
        drr="$(echo "${sra}" | grep -E -o 'DRR[0-9]{4,}' | tr '\n' ',' | sed 's/.$//')"

        biosample="$(esearch -db bioproject -query PRJNA766550 | elink -target biosample | esummary | grep -E -o 'SAMN[0-9]{4,}' | tr '\n' ',' | sed 's/.$//')"

        echo "${assembly};${srr};${err};${drr};${biosample}" || echo "$line failed."
        sleep .5

done 3< rajaAssemblies.txt
