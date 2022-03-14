#!/bin/bash 

# Retrieve all BioSamples from ARGOS BioProject
(esearch -db bioproject -query "PRJEB12890" | elink -target biosample |
	efetch -format xml) > ../PRJEB12890_biosamples.xml