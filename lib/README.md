
## Guide to using scripts in `/lib`:

## biosample_complete.sh

Generates a TSV file (hopefully) adhering to the biosampleMeta schema
as defined in https://github.com/FDA-ARGOS/data.argosdb/blob/main/schema/v1.0/core/biosampleMeta.json

Currently takes a single biosample id argument.  The TSV returned contains
one row per SRA id associated with the biosample.

Example usage:

`./biosample_complete.sh SAMN03255434`

MORE TO COME HERE
