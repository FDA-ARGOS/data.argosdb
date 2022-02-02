#!/bin/bash 

# Retrieve all BioSamples from ARGOS BioProject
(esearch -db bioproject -query "PRJNA231221" | elink -target biosample |
	efetch -format xml) > home/PRJNA231221_biosamples.xml