# HIVE3 Scripts

This directory contains scripts used to generate ARGOS data tables as part of HIVE3-based quality control (QC) workflows.

These scripts were used during earlier ARGOS data releases to process QC outputs, extract relevant metrics, and convert structured JSON outputs into tabular formats for data ingestion.

## Purpose

The scripts in this directory were designed to:
- Process HIVE3 QC output files (e.g., qcAll, ngsQc)
- Extract relevant metrics and metadata
- Convert JSON outputs into TSV-formatted data tables
- Support schema-aligned data generation for ARGOS

**How It Works**
- The pipeline code expects folders of JSON files for each QC category (e.g., qcAll, ngsQc)
- Using the code, these are processed and converted into TSV outputs
- A corresponding columns JSON file (schema-like structure) must be provided for each dataset type in order for the scripts to run correctly

_Each script includes_:
- A header explaining its purpose
- Inline comments describing required inputs and logic
  
## Components
**JSON → TSV Conversion**
- Converts structured HIVE3 output JSON files into tabular TSV format
- Designed to run locally (developed and tested on macOS)
- May require modification for compatibility with Windows environments

**BioSample Data Grabber**
- Retrieves BioSample-associated metadata
- Compatible with updated NCBI data structures as well as HIVE3 outputs
- Refer to script-level documentation for required input criteria

## Important Notes
_Legacy Code_ →
These scripts are no longer actively maintained and may not reflect the current ARGOS pipeline.

_HIVE3 Evolution_ →
HIVE3 outputs and formats changed over time, so some scripts may rely on outdated assumptions.

_Use with Caution_ →
These scripts are best used as references and may require modification to work with current workflows.

## When to Use This Directory
This directory is most useful for:
- Understanding historical ARGOS data processing workflows
- Reviewing how HIVE3 QC outputs were parsed and transformed
- Troubleshooting or comparing against legacy implementations

For actively maintained workflows, see the [../current](lib/current/)
 directory.

## Contact

If you have questions about these scripts, please reach out to the HIVE Lab at mazumder_lab@gwu.edu.
