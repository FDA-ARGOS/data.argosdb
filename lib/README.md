
## Guide to using scripts in `/lib`:

## biosample_complete.sh

Generates a TSV file (hopefully) adhering to the biosampleMeta schema
as defined in https://github.com/FDA-ARGOS/data.argosdb/blob/main/schema/v1.0/core/biosampleMeta.json

Currently requires the `-f` parameter, specifying the path to a text file of bioSample IDs with one ID per line.
The TSV returned contains one row per SRA id associated with the biosample.

Example usage:

`./biosample_complete.sh -f biosample_ids.txt`
