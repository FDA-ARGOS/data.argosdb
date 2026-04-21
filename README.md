# FDA-ARGOS: ARGOSDB

## Introduction
FDA-ARGOS database updates may help researchers rapidly validate diagnostic tests and use qualified genetic sequences to support future product development

As of September 2021, Embleema and George Washington University have been conducting bioinformatic research and system development, focusing on expanding the FDA-ARGOS database. This project expands datasets publicly available in FDA-ARGOS, improves quality control by developing quality matrix tools and scoring approaches that will allow the mining of public sequence databases, and identifies high-quality sequences for upload to the FDA-ARGOS database as regulatory-grade sequences. Building on expansions during the COVID-19 pandemic, this project aims to further improve the utility of the FDA-ARGOS database as a key tool for medical countermeasure development and validation.

For additional details on project information and assembly QC see:

* [ARGOSDB](https://data.argosdb.org/)
* [NCBI Bioproject](https://www.ncbi.nlm.nih.gov/bioproject/231221)
* [Project Information](https://www.fda.gov/emergency-preparedness-and-response/mcm-regulatory-science/expanding-next-generation-sequencing-tools-support-pandemic-preparedness-and-response)

## Data_files
The data files contains raw or processed code. It also contains templates and tables for the data files pushed to data.argosdb.org.
It can be found [here](data_files)
_
## Docs
Documentation can be found in the [docs](docs) directory.

## lib
The [lib](lib/) directory contains scripts and supporting code used for data retrieval, processing, and table generation within ARGOS workflows.

**These scripts have been used to**:
- Retrieve data and identifiers from external APIs (e.g., NCBI)
- Process and transform QC outputs
- Generate standardized data tables and supporting schema inputs

**Subdirectories**

- [HIVE3/](lib/HIVE3/) Contains scripts used historically to generate ARGOS data tables using HIVE3-based QC workflows.
Some scripts may be outdated due to ongoing changes and updates to HIVE3 over time.
- [current/](lib/current/) Contains the actively maintained codebase used to generate the data tables currently available on data.argosdb.org. These scripts can also be adapted to generate custom tables from HIVE QC outputs.
- [shell/](lib/shell/) Contains legacy shell scripts that were used in earlier versions of the pipeline and are no longer actively maintained.

## Schema
This directory contains the ARGOS data schemas, organized by release version. Each version defines the expected structure and validation rules for datasets published on ARGOS.

_Current version_: [v1.6](https://github.com/FDA-ARGOS/data.argosdb/tree/main/schema/v1.6)

**What is included in each schema version?**

Each version contains a core/ subdirectory, which houses JSON schema files corresponding to the data tables available on data.argosdb.org. There is an annotation/ subdirectory, but we are not currently performing any annotations. 

These JSON schemas define:
- Column names (attributes) for each dataset
- Data types (e.g., string, integer, float)
- Descriptions explaining the meaning of each field
- Titles for readability and standardization
- Example values to guide interpretation and usage

This structure ensures consistency, transparency, and reproducibility across ARGOS datasets.

Example schema file:
- [assemblyQC.json](https://github.com/FDA-ARGOS/data.argosdb/blob/main/schema/v1.6/core/assemblyQC.json)
- [ngsQC.json](https://github.com/FDA-ARGOS/data.argosdb/blob/main/schema/v1.6/core/ngsQC.json)

### Validating a data file against a schema:
Assume you wanted to validate a file of the type`SRA_ngsQC`(this same process should work for any of the types we have defined).

- The data file is `/data_files/test_SRA_ngsQC.tsv`
- The schema for a `SRA_ngsQC` data file is `/schema/v0.5/non-core/SRA_ngsQC.json` 

For illitstration purposes cell `T6` in our example data file has been modified. The schema says that the value has to be less than 1, as `gc_ content` is a percentage. The example data sheet has a value of `10.63682374` in that cell, and the following error should be thrown:

`Line 5 failed. '10.63682374' does not match '^[+-]?([0]+\\.?[0-9]*|\\.[0-9]+)$` 

From the project root run:
 
```shell
> python lib/dictionary_utils.py validate -i data_files/test_SRA_ngsQC.tsv -s schema/v0.5/non-core/SRA_ngsQC.json
```

### Validating a data file against a schema with remote files:

Both the schema [`-s`] and input file [`-i`] values can take a URL, assuming they are formatted correctly and resolvable. 

For Example: 
```
>  python lib/dictionary_utils.py validate -i https://raw.githubusercontent.com/FDA-ARGOS/data.argosdb/v0.5/data_files/test_SRA_ngsQC.tsv -s https://raw.githubusercontent.com/FDA-ARGOS/data.argosdb/v0.5/schema/v0.5/non-core/SRA_ngsQC.json
``` 
should give you the same results. 

## covid-drdb
COVID-DRDB is created by the HIVDB team of Stanford University. It includes resistance data of SARS-CoV-2 for convalescent plasma, vaccinee plasma and monoclonal antibodies collected from published peer-reviewed/pre-print studies.
The documents and files can be found [here](covid-drdb/README.md).

## Citations
Sichtig, H., Minogue, T., Yan, Y. et al. FDA-ARGOS is a database with public quality-controlled reference genomes for diagnostic use and regulatory science. Nat Commun 10, 3313 (2019). https://doi.org/10.1038/s41467-019-11306-6


