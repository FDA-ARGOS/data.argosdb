# lib

The lib directory contains scripts and supporting code used for data retrieval, processing, validation, and table generation within ARGOS workflows.

These scripts have been used to:
- Retrieve data and identifiers from external APIs (e.g., NCBI)
- Process and transform QC outputs
- Generate standardized data tables
- Validate datasets against defined schemas

⚠️ Note: Much of the code in this directory reflects earlier versions of the pipeline and may no longer be part of the current workflow. Except for the current or HIVE3 olders. 

## Subdirectories
`HIVE3/`
Scripts used historically to generate ARGOS data tables using HIVE3-based QC workflows. Some scripts may be outdated due to changes in HIVE3 over time.

`current/`
Contains the actively maintained codebase used to generate the data tables currently available on data.argosdb.org.

`shell/`
Legacy shell scripts used in earlier versions of the pipeline. These are no longer actively maintained.

## Notes
Many scripts in this directory are deprecated or legacy

Retained for:
- Reproducibility
- Historical reference
- Troubleshooting older datasets

For current workflows, refer to the main README and updated pipelines.

## Biosample Schema script
A much older file that is still applicable, biosample_complete.sh, contains the steps below.

Generates a TSV file intended to follow the biosampleMeta schema:
https://github.com/FDA-ARGOS/data.argosdb/blob/main/schema/v1.4/core/biosampleMeta_HIVE.json

Output: One row per SRA ID associated with each biosample
Required parameters:
`-f`: Path to text file of BioSample IDs (one per line)
`-n`: Path to NGS QC file (TSV format)

Optional parameters:
`-b`: BCO ID
`-s`: Schema version
`-d`: Debug mode

Example:
`./biosample_complete.sh -f biosample_ids.txt -n ngsQC_HIVE.tsv -b ARGOS_000028 -s v1.12`

## Schema Validation
The dictionary_utils.py script can be used to validate data tables against ARGOS schemas. This is an older file, but can still be applicable.

_Local validation_
`python lib/dictionary_utils.py validate \
  -i data_files/test_SRA_ngsQC.tsv \
  -s schema/v0.5/non-core/SRA_ngsQC.json`

_Remote validation_
`python lib/dictionary_utils.py validate \
  -i https://raw.githubusercontent.com/FDA-ARGOS/data.argosdb/v0.5/data_files/test_SRA_ngsQC.tsv \
  -s https://raw.githubusercontent.com/FDA-ARGOS/data.argosdb/v0.5/schema/v0.5/non-core/SRA_ngsQC.json`

_Example error_
`Line 5 failed. '10.63682374' does not match '^[+-]?([0]+\\.?[0-9]*|\\.[0-9]+)$'`
