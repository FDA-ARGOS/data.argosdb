#!/bin/bash 

# # Retrieve all run IDs from an ARGOS BioSample
esearch -db biosample -query $1 | elink -target sra | efetch -format docsum
