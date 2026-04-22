# Current Data Table Scripts
This directory contains the actively maintained scripts used to generate the data tables available on data.argosdb.org.

These scripts represent the current ARGOS workflow for transforming QC outputs and associated metadata into standardized, schema-aligned data tables.

## Overview
The scripts in this directory are used to:
- Process QC outputs (e.g., NGS QC data)
- Retrieve and integrate BioSample metadata
- Generate finalized TSV ngs, assembly, and biosample data tables for ARGOS
**Last updated**: May 22, 2025

## Recommended Scripts
To ensure compatibility with the current data structure and workflows, use the following:

- `biosample_metadata_grabberV4`
   - Use this version for all current table generation
   - Updated to handle the removal of assembly IDs from JSON inputs (previous JSONs did not have this)
   - Replaces previous versions of the BioSample metadata pipeline
- `ngsQC_datatable_V2.sh`
   - Use this script for generating NGS QC data tables
   - Includes updated computations and formatting logic
   - Earlier versions are deprecated and may not function correctly

## Versioning Notes
Previous versions of these scripts (e.g., biosample_metadata_grabber_V3, older NGS QC scripts) have been moved to the [../lib](lib/) or [../HIVE3](lib/HIVE3/) directory for reference. Only the scripts in this directory are considered _current_ and supported. 

## Usage Notes
Scripts are designed to work with ARGOS QC outputs and expected input formats. Ensure that input data aligns with the latest schema definitions. Refer to inline script documentation for required inputs and parameters.
