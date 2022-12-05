#!/bin/bash

# Curl for NCBI Trace to retrieve SRR analysis values
# Need to provide a text file with a list of SRR accessions
#
# $1 is path to text file with list of SRR accessions
# $2 is path to folder to create with SRR analysis XMLs

mkdir -p $2

for i in $(cat $1)
do
         echo "Fetching $i"
         (esearch -db sra -query $i | efetch -format xml) > $2/$i.xml
done
