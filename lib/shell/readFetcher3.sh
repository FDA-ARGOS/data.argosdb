while read -u3 row
do

    echo "Reading $row..."

	assembly="$(esearch -db sra -query $row | elink -target assembly | esummary)" | grep "FtpPath_RefSeq" | sed -r 's/.*>(ftp:\/\/ftp.ncbi.nlm.nih.gov\/genomes\/all\/.+.+)<.*/\1/' | rev | cut -d'/' -f 1 | rev >> assemblies.txt || echo "$row failed."
    sleep .4

done 3< readList3.txt