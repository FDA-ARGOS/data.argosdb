#!/bin/bash
# Pull all assemblies associated with a list of BioProjects, preferentially pulling RefSeq assemblies if they exist, then falling back to GenBank, if not.
# This script is written for anonymous access (sleep .8 seconds per query) since it's a public repo (and I haven't set up .secrets yet). Use your NCBI API auth key, if you like.

while read -u3 line
do
    if [ -z "$(esearch -db bioproject -query $line | elink -target assembly | esummary | grep "FtpPath_RefSeq" | sed -r 's/.*>(ftp:\/\/ftp.ncbi.nlm.nih.gov\/genomes\/all\/.+.+)<.*/\1/' | rev | cut -d'/' -f 1 | rev >>assemblies.txt ; sleep .8)" ]; then
        echo "No RefSeq entries found for $line, trying GenBank..."
        esearch -db bioproject -query $line | elink -target assembly | esummary | grep -i 'FtpPath type="GenBank"' | sed -r 's/.*>(ftp:\/\/ftp.ncbi.nlm.nih.gov\/genomes\/all\/.+.+)<.*/\1/' | rev | cut -d'/' -f 1 | rev >>assemblies.txt ; sleep .8
    else
        echo "RefSeq entries found for $line"
    fi
done 3< rajaAssemblies.txt
