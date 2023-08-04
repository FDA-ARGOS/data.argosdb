#!/bin/bash 
#
# Usage: ./biosample_complete.sh -f BIOSAMPLEIDFILE -n NGSQCFILE -b BCO_ID -s SCHEMA_VERSION


WORKING_DIR=biosample_$(date "+%Y%m%d_%H%M%S")
mkdir -p $WORKING_DIR

while getopts "f:n:b:s:" opt; do
  case $opt in
    f | file)
      sampleidfile=$OPTARG
      ;;
    n | ngsqc_file )
      ngsqcfile=$OPTARG
      ;;
    b | bco_id )
      bco_id=$OPTARG
      ;;
    s | schema )
      schema=$OPTARG
      ;;
  esac
done

# Read sample id file into an array
IDLINES=($(cat $sampleidfile))

for biosample_id in ${IDLINES[@]}
do
  echo "Retrieving biosample metadata for $biosample_id"
  # Retrieve biosample metadata
  esearch -db biosample -query $biosample_id | efetch -format xml > $WORKING_DIR/biosample_meta_$biosample_id.xml
  # parse metadata
  echo ...Parsing XML biosample metadata to create tsv
  python3 ./sra_sample_parser.py -n $ngsqcfile -f $WORKING_DIR/biosample_meta_$biosample_id.xml -b $bco_id -s $schema > $WORKING_DIR/biosample_meta_$biosample_id.tsv
done 

echo Biosample metadata files are located at $WORKING_DIR
 
TSVFILENAME=$WORKING_DIR/${WORKING_DIR}_meta_combined.tsv

echo Adding header line
head -n 1 $WORKING_DIR/biosample_meta_${IDLINES[0]}.tsv > $TSVFILENAME

echo looping through lines
for biosample_id in ${IDLINES[@]}
do
   tail -n +2 $WORKING_DIR/biosample_meta_$biosample_id.tsv >> $TSVFILENAME
done

echo Combined tsv file is located at $TSVFILENAME
