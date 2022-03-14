#!/bin/bash 

# # Retrieve all run IDs from ARGOS BioProject
(esearch -db bioproject -query "PRJEB12890" | elink -target sra | \
	efetch -format docsum) > home/PRJEB12890_runs.xml