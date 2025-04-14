#Christie Woodside
#April 14 2025

'''Counts the number of SRR IDs in each of the kingdoms from the ARGOS bioProject specifically. Takes in the SRA total entires summary file. 
To download click on the SRR hyperlink from the bioproject > send to > file > summary > download. The count will be displayed in the terminal'''


import csv
import argparse
import os
import time
import requests
import xml.etree.ElementTree as ET
from collections import defaultdict


def get_kingdom_from_organism(organism):
    api_key = os.getenv("NCBI_API_KEY") or "your API key from NCBI"            # <---------- make sure to add your NCBI API key here

    # Step 1: Get Taxonomy ID
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_params = {
        "db": "taxonomy",
        "term": organism,
        "retmode": "json",
        "api_key": api_key
    }
    search_resp = requests.get(search_url, params=search_params)
    if search_resp.status_code == 429:
        time.sleep(1)
        search_resp = requests.get(search_url, params=search_params)
    if search_resp.status_code != 200:
        return None

    try:
        tax_id = search_resp.json()["esearchresult"]["idlist"][0]
    except Exception:
        return None

    # Step 2: Get lineage
    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    fetch_params = {
        "db": "taxonomy",
        "id": tax_id,
        "retmode": "xml",
        "api_key": api_key
    }
    fetch_resp = requests.get(fetch_url, params=fetch_params)
    if fetch_resp.status_code == 429:
        time.sleep(1)
        fetch_resp = requests.get(fetch_url, params=fetch_params)
    if fetch_resp.status_code != 200:
        return None

    try:
        root = ET.fromstring(fetch_resp.content)
        lineage = root.findtext("Taxon/Lineage")
        if "Viruses" in lineage:
            return "virus"
        elif "Fungi" in lineage:
            return "fungi"
        elif "Bacteria" in lineage:
            return "bacteria"
        else:
            return "other"
    except Exception:
        return None

def main():
    parser = argparse.ArgumentParser(description="Count SRR entries by kingdom from a CSV file")
    parser.add_argument("input_file", help="Path to the CSV file")
    args = parser.parse_args()

    kingdom_counts = defaultdict(int)
    organism_to_kingdom = {}
    skipped_blank = 0
    failed_taxonomy = 0
    counted_total = 0
    other_kingdom = 0

    with open(args.input_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        entries = list(reader)

        for i, row in enumerate(entries, 1):
            organism = row.get("Organism Name", "").strip()
            if not organism:
                print("------------empty: ", i)
                skipped_blank += 1
                continue

            if organism not in organism_to_kingdom:
                kingdom = get_kingdom_from_organism(organism)
                if kingdom is None:
                    print("Organism failed taxonomy: ", organism)
                    failed_taxonomy += 1
                    continue
                organism_to_kingdom[organism] = kingdom
            else:
                kingdom = organism_to_kingdom[organism]

            if kingdom in {"fungi", "bacteria", "virus"}:
                kingdom_counts[kingdom] += 1
                counted_total += 1
            else:
                print("Organism untracked: ", organism)
                other_kingdom += 1

            #print(f"[{i}/{len(entries)}] {organism} → {kingdom}")

    print("\nTotal SRR counts by kingdom:")
    print(f"fungi\t{kingdom_counts['fungi']}")
    print(f"bacteria\t{kingdom_counts['bacteria']}")
    print(f"virus\t{kingdom_counts['virus']}")

    print(f"\nSummary:")
    print(f"Total rows in CSV: {len(entries)}")
    print(f"Counted entries: {counted_total}")
    print(f"Skipped (blank organism): {skipped_blank}")
    print(f"Failed taxonomy lookup: {failed_taxonomy}")
    print(f"Other/untracked kingdom: {other_kingdom}")

if __name__ == "__main__":
    main()
