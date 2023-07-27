
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

Example usage:

`./biosample_complete.sh -f biosample_ids.txt -n ngsQC_HIVE.tsv -b ARGOS_000028 -s v1.12`
