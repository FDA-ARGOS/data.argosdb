# data.argosdb

Repository for data files and schema definitions for ARGOS rpoject. 

## Data_files
Space for raw or processed data files
_
## Docs
space for documentation

## Examples
Space for example data files

## Schema
Houses the ARGOS data schemas

## lib
For scripts and sutch


### Too validate against a Schema

```shell
python lib/validate.py 

usage: argosdb [options]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -j JSON, --json JSON  JSON to process.
  -s SCHEMA, --schema SCHEMA
                        Root json schema to validate against.

```