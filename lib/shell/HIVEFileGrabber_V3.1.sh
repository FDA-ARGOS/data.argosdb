#!/usr/bin/env bash

#Remove previous output, if it exists, and create a new one:
FILE=/home/transfer_hive_prod/outputList.txt
if test -f "$FILE"; then
	echo "previous output file exists, removing."
	rm /home/transfer_hive_prod/outputList.txt || echo "Can't rm output."
fi
touch /home/transfer_hive_prod/outputList.txt || echo "Can't touch output."


#Set directories and IO:
ROOTDIR=~/data/temp
INFILE=$ROOTDIR/inputList.txt
OUTFILE=/home/transfer_hive_prod/outputList.txt

#Remove commas from list of input HIVE IDs:
sed -i 's/,//g' $INFILE || echo "Can't sed 1"

#Remove list numbers/colons from list of input HIVE IDs:
sed -i 's/[[:digit:]]\{1,\}:[[:space:]]//g' $INFILE || echo "Can't sed 2"

#Print column headers:
echo "file_name;codon_table;%protein_coding;%not_coding;%_reads;density_Ns_per_read;%complex;%not_complex;avg_quality_a;avg_quality_t;avg_quality_g;avg_quality_c;count_a;count_c;count_g;count_t" >> /home/transfer_hive_prod/outputList.txt || echo "Can't set headers."

#For every input object ID, go to the appropriate directory in the back end of HIVE2, get the complexity, codonQC, sumLetter, and countNs tables and append their contents to a merged file:
cat $INFILE | while read line
do
	echo "Processing $line..."
	qcd $line || echo "Can't qcd."

#Test to see if the line is an SRR or an ERR:
	prefixTest=$(grep --text 'SRR.*fastq' ./_.zip)
	exitStatus=$?
	if [[ $exitStatus == 0 ]]
        then
		grep -o --text 'SRR.*fastq' ./_.zip | sed '1d' > ./TEMPIDFile.txt || echo "Can't grep SRR."
	else
		echo "No SRR pattern detected, trying ERR..."
		grep -o --text 'ERR.*fastq' ./_.zip | sed '1d' > ./TEMPIDFile.txt || echo "Can't grep ERR."
	fi

#Take only the data (minus the header) and store it in a temporary file:
	cat _.qc2.codonQCTable.csv | sed '1d' > ./TEMPcodonQCTable.csv || echo "Can't cat codonQCTable"
	cat _.qc2.sumLetterTable.csv > ./TEMPLetterCountQuality.csv || echo "Can't cat sumLetterTable"
	cat _.qc2.countNsPercentageTable.csv | sed '1d' > ./TEMPcountNsPercentageTable.csv || echo "Can't cat countNsTable"
	cat _.qc2.ComplexityTable.csv > ./TEMPComplexityTable.csv || echo "Can't cat ComplexityTable"

#Combine and format all of the output (consider just making this an awk in the future):
	python3 $ROOTDIR/QCTableFormatter3.1.py >> /home/transfer_hive_prod/outputList.txt || echo "Can't python"

#Clean up:
	rm ./TEMPIDFile.txt || echo "Can't rm file 1"
	rm ./TEMPcodonQCTable.csv || echo "Can't rm file 2"
	rm ./TEMPLetterCountQuality.csv || echo "Can't rm file 3"
	rm ./TEMPcountNsPercentageTable.csv || echo "Can't rm file 4"
	rm ./TEMPComplexityTable.csv || echo "Can't rm file 5"

done

echo "Done."
