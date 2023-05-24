#!/bin/bash 

# # Retrieve all run IDs from ARGOS BioProject
(esearch -db bioproject -query "PRJNA231221" | elink -target sra | \
	efetch -format docsum) > home/ngsqc/PRJNA231221_runs.xml