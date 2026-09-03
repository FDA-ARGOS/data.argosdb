#!/usr/bin/env bash

# Identify SRA run accessions and linked RefSeq assemblies for a list
# of NCBI BioProject accessions.
#
# Usage:
#   ./find_linked_sra_assemblies.sh \
#       bioprojects.txt \
#       read_accessions.txt \
#       assembly_accessions.txt
#
# BioSample records are not retrieved because they are collected later
# by the HIVE workflow.
#
# Requirements:
#   NCBI Entrez Direct: esearch, elink, esummary, and xtract
#
# NCBI_API_KEY may be set in the environment for authenticated requests.

set -uo pipefail

readonly INPUT_FILE="${1:-rajaBioprojects.txt}"
readonly RUN_OUTPUT="${2:-read_accessions.txt}"
readonly ASSEMBLY_OUTPUT="${3:-assembly_accessions.txt}"
readonly REQUEST_DELAY="${NCBI_REQUEST_DELAY:-0.4}"

for required_command in esearch elink esummary xtract sort mktemp; do
    if ! command -v "${required_command}" >/dev/null 2>&1; then
        echo "Error: required command not found: ${required_command}" >&2
        exit 127
    fi
done

if [[ ! -r "${INPUT_FILE}" ]]; then
    echo "Error: input file is not readable: ${INPUT_FILE}" >&2
    exit 1
fi

temporary_directory="$(mktemp -d)"
trap 'rm -rf "${temporary_directory}"' EXIT

raw_runs="${temporary_directory}/runs.txt"
raw_assemblies="${temporary_directory}/assemblies.txt"

: > "${raw_runs}"
: > "${raw_assemblies}"

while IFS= read -r bioproject || [[ -n "${bioproject}" ]]; do
    # Ignore blank lines and comments.
    [[ -z "${bioproject}" || "${bioproject}" == \#* ]] && continue

    echo "Finding runs associated with ${bioproject}..." >&2

    if ! sra_summary="$(
        esearch -db bioproject -query "${bioproject}" |
            elink -target sra |
            esummary
    )"; then
        echo "Warning: NCBI query failed for ${bioproject}." >&2
        continue
    fi

    run_accessions="$(
        printf '%s\n' "${sra_summary}" |
            grep -Eo '(SRR|ERR|DRR)[0-9]+' |
            sort -u ||
            true
    )"

    if [[ -z "${run_accessions}" ]]; then
        echo "Warning: no SRA runs found for ${bioproject}." >&2
    else
        printf '%s\n' "${run_accessions}" >> "${raw_runs}"
    fi

    sleep "${REQUEST_DELAY}"
done < "${INPUT_FILE}"

# Overwrite rather than append so that rerunning the script produces the
# same output for the same NCBI results.
sort -u "${raw_runs}" > "${RUN_OUTPUT}"

while IFS= read -r run_accession || [[ -n "${run_accession}" ]]; do
    [[ -z "${run_accession}" ]] && continue

    echo "Finding assemblies associated with ${run_accession}..." >&2

    if ! assembly_summary="$(
        esearch -db sra -query "${run_accession}" |
            elink -target assembly |
            esummary
    )"; then
        echo "Warning: NCBI query failed for ${run_accession}." >&2
        continue
    fi

    assembly_accessions="$(
        printf '%s\n' "${assembly_summary}" |
            xtract -pattern DocumentSummary -element AssemblyAccession |
            awk 'NF' |
            sort -u ||
            true
    )"

    if [[ -z "${assembly_accessions}" ]]; then
        echo "Warning: no linked assembly found for ${run_accession}." >&2
    else
        printf '%s\n' "${assembly_accessions}" >> "${raw_assemblies}"
    fi

    sleep "${REQUEST_DELAY}"
done < "${RUN_OUTPUT}"

sort -u "${raw_assemblies}" > "${ASSEMBLY_OUTPUT}"

echo "Wrote $(wc -l < "${RUN_OUTPUT}") unique run accessions to ${RUN_OUTPUT}." >&2
echo "Wrote $(wc -l < "${ASSEMBLY_OUTPUT}") unique assembly accessions to ${ASSEMBLY_OUTPUT}." >&2
