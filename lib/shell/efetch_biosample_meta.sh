#!/bin/bash 

# # Retrieve biosample metadata for a single biosample accession number
esearch -db biosample -query $1 | efetch -format xml
