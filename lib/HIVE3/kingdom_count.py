# Christie Woodside
# April 14, 2025

'''Counts the number of biosamples associated with each kingdom for the FDA-ARGOS BioProject. This will be used for the table in the paper. Takes in the biosample summary text file.
To doanload, click on the biosample hyperlink in the bioproject page > send to > file > summary(text) > download. It should be a .txt file'''

import re
import requests
import argparse
import os
import time
import xml.etree.ElementTree as ET
from collections import Counter

def parse_biosample_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    entries = re.split(r'\n(?=\d+\. Pathogen:)', content)
    organisms = []

    for entry in entries:
        organism_match = re.search(r'Organism:\s*(.+)', entry)
        if organism_match:
            organisms.append(organism_match.group(1).strip())

    return organisms

def get_kingdom(organism):
    api_key = os.getenv("NCBI_API_KEY") or "your own API key"                    #<------------ make sure to add your NCBI API key to this

    # Get tax ID
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_params = {
        "db": "taxonomy",
        "term": organism,
        "retmode": "json",
        "api_key": api_key
    }
    search_resp = requests.get(esearch_url, params=search_params)
    if search_resp.status_code == 429:
        time.sleep(1)
        search_resp = requests.get(esearch_url, params=search_params)
    if search_resp.status_code != 200:
        return "Unknown"

    try:
        tax_id = search_resp.json()["esearchresult"]["idlist"][0]
    except Exception:
        return "Unknown"

    # fetch the lineage
    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    fetch_params = {
        "db": "taxonomy",
        "id": tax_id,
        "retmode": "xml",
        "api_key": api_key
    }
    fetch_resp = requests.get(efetch_url, params=fetch_params)
    if fetch_resp.status_code == 429:
        time.sleep(1)
        fetch_resp = requests.get(efetch_url, params=fetch_params)
    if fetch_resp.status_code != 200:
        return "Unknown"

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
            return "Other"
    except Exception:
        return "Unknown"

def main():
    parser = argparse.ArgumentParser(description="Count BioSamples by Kingdom")
    parser.add_argument("input_file", help="Path to biosample_result.txt")
    args = parser.parse_args()

    organisms = parse_biosample_file(args.input_file)

    kingdom_counts = Counter()

    for i, organism in enumerate(organisms, 1):
        kingdom = get_kingdom(organism)
        kingdom_counts[kingdom] += 1
        #print(f"[{i}/{len(organisms)}] {organism} => {kingdom}")

    print("\nBiosample counts by kingdom:")
    print(f"fungi\t{kingdom_counts.get('fungi', 0)}")
    print(f"bacteria\t{kingdom_counts.get('bacteria', 0)}")
    print(f"virus\t{kingdom_counts.get('virus', 0)}")

if __name__ == "__main__":
    main()
