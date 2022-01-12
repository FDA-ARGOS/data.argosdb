#!/bin/bash

# Curl for NCBI Trace to retrieve SRR analysis values
# Need to provide a text file with a list of SRR accessions
for i in $(cat home/PRJNA231221_runs.txt )
do
	 var1='https://trace.ncbi.nlm.nih.gov/Traces/sra/?run='
	 var2=$i'&retmode=xml'
	 var3=$var1$var2
	 curl -o home/SRR_xmls/$i.xml $var3
done
