#!/usr/bin/env bash

#BioSample metadata grabber. This script pull data from an NCBI biosample record and formats it into ARGOS QC table format.
#This script assumes that eutils is already installed for the uers.

#Get accession number from user:
echo ""
echo "------------------------------------------------------------------------------------------------------------"
echo "| This script is limited to one accession at a time to prevent banning from NCBI servers for rapid queries.|"
echo "| Results will be stored in /home/keeneyjg/projects/ARGOS/data/biosample_record.csv                        |"
echo "------------------------------------------------------------------------------------------------------------"
echo ""
read -p "Enter accession number (in the form of SAMN########):" ACCESSION || echo "Can't read input."

#Check for proper input:
if [[ $ACCESSION =~ ^SAMN[0-9]{8}$ ]];
then

#Remove previous outfile if it exists and set output locations.
rm /home/keeneyjg/projects/ARGOS/data/biosample_record.csv || echo "Didn't remove outfile."
OUTFILE="/home/keeneyjg/projects/ARGOS/data/biosample_record.csv"
TEMPFILE="/home/keeneyjg/projects/ARGOS/data/tempFile.txt"
touch $OUTFILE || echo "Can't touch outfile."
touch $TEMPFILE || echo "Can't touch tempfile."

#Storing this command here in case needed; this is the EBI version via their API:
#curl 'https://www.ebi.ac.uk/biosamples/samples/SAMEA3219152/curationlinks' -i -X GET -H 'Content-Type: application/json'
#For more, see: https://www.ebi.ac.uk/biosamples/docs/

#Use eutils to grab the biosample record based on user input:
efetch -db biosample -id $ACCESSION > $TEMPFILE || echo "Can't eutils."

#Collect relevant fields from biosample record, set equal to - if no data:
NAME=$(sed -n -e 's/^Organism: //p' $TEMPFILE)
if [ -z "$NAME" ]
then
$NAME="-"
fi

STRAIN=$(sed -n -e 's/[[:space:]]*\/strain="//p' $TEMPFILE)
if [ -z "$STRAIN" ]
then
$STRAIN="-"
fi


#Print collected output:
echo "Header1;Header2" >> $OUTFILE || echo "Can't set headers."
echo "$NAME;$STRAIN" >> $OUTFILE || echo "WARNING: didn't print all fields."

#Clean up.
rm $TEMPFILE || echo "Can't remove tempfile."

else
echo "Wrong input format detected, exiting."
exit
fi

