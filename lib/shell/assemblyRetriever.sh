#!/bin/bash
echo "Assembly;SRR;ERR;DRR;BioSample" > ~/ukAssemblies.txt

# Get all assemblies for each BioProject:
while read -u3 row
do

        assembly="$(esearch -db bioproject -query $row | elink -target assembly | esummary | grep "FtpPath_RefSeq" | sed -r 's/.*>(ftp:\/\/ftp.ncbi.nlm.nih.gov\/genomes\/all\/.+.+)<.*/\1/' | rev | cut -d'/' -f 1 | rev | tr '\n' ',' | sed 's/.$//')" >> ~/assemblies.tmp

	echo "${assembly}" >> ~/assemblies.tmp

	sleep .5
done 3< ~/workspace/ukBioprojects/rajaAssemblies.txt

# Get all read data and BioSample information for each assembly:
while read -u3 line
do
	echo "Using $line"

        sra="$(esearch -db bioproject -query $line | elink -target sra | esummary)"
	echo "Using "${sra}"
        srr="$(echo "${sra}" | grep -E -o 'SRR[0-9]{4,}' | tr '\n' ',' | sed 's/.$//')"
        err="$(echo "${sra}" | grep -E -o 'ERR[0-9]{4,}' | tr '\n' ',' | sed 's/.$//')"
        drr="$(echo "${sra}" | grep -E -o 'DRR[0-9]{4,}' | tr '\n' ',' | sed 's/.$//')"

        #biosample="$(esearch -db bioproject -query $line | elink -target biosample | esummary | grep -E -o 'SAMN[0-9]{4,}' | tr '\n' ',' | sed 's/.$//')"

        echo "$line;${srr};${err};${drr}" || echo "$line failed." >> ~/ukAssemblies.txt
        sleep .5

done 3< ~/assemblies.tmp
