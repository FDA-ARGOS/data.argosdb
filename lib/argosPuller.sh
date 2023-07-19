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
