#!/bin/bash 

WORKING_DIR=biosample_$(date "+%Y%m%d_%H%M%S")
mkdir -p $WORKING_DIR

# # Retrieve biosample metadata
echo Retrieving biosample metadata for $1
./efetch_biosample_meta.sh $1 > $WORKING_DIR/biosample_meta_$1.xml

# parse metadata
echo Parsing XML biosample metadata to create tsv
python3 ./sra_sample_parser.py -f $WORKING_DIR/biosample_meta_$1.xml > $WORKING_DIR/biosample_meta_$1.tsv

# # Retrieve all run IDs from an ARGOS BioSample
echo Retrieving run data associated with biosample $1
./efetch_biosample_runs.sh $1 > $WORKING_DIR/biosample_runs_$1.xml

# Get experiment IDs list
echo Parsing run data to extract run IDs
python3 ./xml_feature_parser.py -f $WORKING_DIR/biosample_runs_$1.xml -o $WORKING_DIR/biosample_runlist_$1.txt

# Fetch experiment info based on run list, place in $WORKING_DIR/biosample_runs_$1
echo Fetching experiment metadata assoicated with each run ID
./eutils_SRAfetch_experiments.sh $WORKING_DIR/biosample_runlist_$1.txt $WORKING_DIR/biosample_runs_$1 
# Parse experiment info
echo Parsing experiment metadata XMLs to create tsv with metadata from runs
python3 ./sra_run_parser.py -d $WORKING_DIR/biosample_runs_$1/ > $WORKING_DIR/biosample_run_$1.tsv

# Merge TSV files
echo Merging biosample metadata and run metadata tsv files
python3 ./biosample_tsv_join.py -m $WORKING_DIR/biosample_meta_$1.tsv -r $WORKING_DIR/biosample_run_$1.tsv -o $WORKING_DIR/biosample_complete_$1.tsv

echo Biosample metadata tsv file is located at $WORKING_DIR/biosample_complete_$1.tsv
