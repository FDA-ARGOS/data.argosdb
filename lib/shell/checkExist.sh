#!/usr/bin/env bash

# With large batches, it was discovered that some quitely fail. Cause is unknown.
# This script just iterates through all objectIDs (provided in a separate file) to find the ones that quietly failed.
# Run on HIVE3 as admin.

INFILE=objectIDNewline.txt

cat $INFILE | while read line
do
#	echo "Processing $line..."
	qcd $line || echo "Can't qcd $line."
	if [ -e -qcAll.json ]
	then
		echo "$line OK."
	else
		echo "~~~~~~~~CANNOT FIND $line!~~~~~~~~~"
	fi
done
