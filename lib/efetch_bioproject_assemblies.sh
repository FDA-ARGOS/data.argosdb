#!/bin/bash 

# # Retrieve all assemblies from ARGOS BioProject
(esearch -db bioproject -query "PRJNA231221" | elink -target assembly |
	efetch -format docsum) > PRJNA231221_assemblies.xml