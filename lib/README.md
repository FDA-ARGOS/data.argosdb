
## Guide to using scripts in `/lib`:

## biosample_complete.sh

Generates a TSV file (hopefully) adhering to the biosampleMeta schema
as defined in https://github.com/FDA-ARGOS/data.argosdb/blob/main/schema/v1.4/core/biosampleMeta_HIVE.json

The TSV returned contains one row per SRA id associated with the biosample.

Required parameters:
- `-f`: Path to text file of bioSample IDs (one ID per line)
- `-n`: Path to NGS QC file (must be TSV format)

Optional parameters:
- `-b`: BCO ID
- `-s`: schema version
- `-d`: debug (if `-d` is not set, the output directory will contain intermediate xml and tsv files)

Example usage:

`./biosample_complete.sh -f biosample_ids.txt -n ngsQC_HIVE.tsv -b ARGOS_000028 -s v1.12`

## lib
For scripts and sutch


### Validating a data file against a schema:
Assume you wanted to validate a flie of the type`SRA_ngsQC`(this same process should work for any of the types we have defined).

- The data file is `/data_files/test_SRA_ngsQC.tsv`
- The schema for a `SRA_ngsQC` data file is `/schema/v0.5/non-core/SRA_ngsQC.json` 

For illitstration purposes cell `T6` in our example data file has been modified. The schema says that the value has to be less than 1, as `gc_ content` is a percentage. The example data sheet has a value of `10.63682374` in that cell, and the following error shoudl be thrown:

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

