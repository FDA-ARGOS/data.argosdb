#!/bin/bash
# Script to look through our processed data to find those requested by BARDA. This file looks through the assemblyQC_ARGOS" file.
while read line; do
        echo "Total number of entries for: $line:"
        cat /data/shared/argosdb/generated/datasets/reviewed/assemblyQC_ARGOS_extended.tsv | awk -F '\t' -v name="$line" '$1 ~ name' | wc -l
        echo -e ' \t '"Unique assembles:"
        echo -ne ' \t ' && cat /data/shared/argosdb/generated/datasets/reviewed/assemblyQC_ARGOS_extended.tsv |  awk -F '\t' -v name="$line" '$1 ~ name' | awk -F '\t' '{ print $4 }' | sort | uniq | wc -l
done < bardaList.txt
