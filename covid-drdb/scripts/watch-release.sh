#! /bin/bash

set -e

echo 'Waiting for the release workflow to start...'
sleep 15

gh run watch --repo hivdb/covid-drdb-payload $(gh run list --repo hivdb/covid-drdb-payload --json headBranch,databaseId,workflowName -L 4 --jq ".[] | select( .headBranch == \"$1\" and .workflowName == \"Release & pre-release\").databaseId")
