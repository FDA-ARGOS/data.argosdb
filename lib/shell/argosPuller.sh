#!/bin/bash

if [ "$#" -eq 0 ]
then
        echo "No arguments supplied"
        exit 1
fi

for FILE in /software/argosdb/dataset-maker/recipes/*.json;
do
        query=$(printf "%s\n" "${FILE##*/}" | sed -r 's/(.*[0-9]{6})_.*/\1/' | tr '[:lower:]' '[:upper:]' $sample);
        echo $query;
        if [[ ! -e ~/api/data/$query.json ]]; then
                touch ~/api/data/"$query".json;
        fi;
        curl -X 'POST' \
                'https://api.argosdb.org/dataset/detail' \
                -H 'accept: application/json' \
                -H 'Content-Type: application/json' \
                -d '{
                    "bcoid": "'"$query"'",
                    "dataversion": "'"$1"'"
                }' > ~/api/data/"$query".json;
done

echo "Finished downlaoding data from ARGOS API. Attempting to push data sets to HIVE 3...\n"
read -p "Enter your login for HIVE 3: " user
read -sp "Enter your password" pass
curl -k -c ~/hive3.cookie "https://hive3.biochemistry.gwu.edu/dna.cgi?cmdr=login&login=${user}&pswd=${pass}"

for DATASET in ~/api/data/*.json;
#for DATASET in ~/api/data/*.json;
do
        curl -v -b  ~/hive3.cookie  -X POST -F "content=@${DATASET}" -F 'cmd=objSetFile' -F 'raw=1'  -F 'bin=0' -F 'type=u-file' -F 'filename=_.json' -F 'name=apiTestFile'  -F 'ext=json'  'https://hive3.biochemistry.gwu.edu/dna.cgi';
done




