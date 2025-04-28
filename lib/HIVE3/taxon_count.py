#Christie Woodside
#April 14 2025

'''Code will count the number of unique taxon IDs from the bioproject so we get an accurate 
count of bacteria, fungi, and viruses that are therefore in the bioproject. Takes in the biosample file containing the biosample summaries from the bioproject. To download
click on the biosample hyperlink from the bioproject > send to > file > summary(text) > download'''

# python3 taxon_count.py /Users/user/Desktop/biosample_result.txt
import re
import requests
import argparse
import os
import time
import xml.etree.ElementTree as ET
from collections import defaultdict

# def parse_biosample_file(filepath):
#     with open(filepath, 'r') as f:
#         content = f.read()

#     entries = re.split(r'\n(?=\d+\. Pathogen:)', content)
#     organisms = []

#     for entry in entries:
#         organism_match = re.search(r'Organism:\s*(.+)', entry)
#         if organism_match:
#             organisms.append(organism_match.group(1).strip())

#     return organisms

def parse_biosample_file(filepath):
    organisms = []

    with open(filepath, 'r') as f:
        for idx, line in enumerate(f, 1):
            line = line.strip()
            if line.startswith("Organism:"):
                organism = line.replace("Organism:", "").strip()
                if organism:
                    organisms.append(organism)
                else:
                    print(f"⚠️ Warning: Empty organism at line {idx}")

    print(f"Parsed {len(organisms)} organisms from file.")
    return organisms

def get_tax_id_and_kingdom(organism):
    api_key = os.getenv("NCBI_API_KEY") or "YOUR API KEY"    #<----- your API key goes here

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
        print('--- bad status:', organism)
        return None, None

    try:
        tax_id = search_resp.json()["esearchresult"]["idlist"][0]
    except Exception:
        print("-- no id list found for: ", organism)
        return None, None

    # Step 2: Get kingdom from lineage
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
        print('--- bad status code again: ', organism)
        return tax_id, None

    try:
        root = ET.fromstring(fetch_resp.content)
        lineage = root.findtext("Taxon/Lineage")
        if "Viruses" in lineage:
            return tax_id, "virus"
        elif "Fungi" in lineage:
            return tax_id, "fungi"
        elif "Bacteria" in lineage:
            return tax_id, "bacteria"
        else:
            print('--- for Other org:   ', lineage)
            return tax_id, "other"
    except Exception:
        return tax_id, "Unknown"

def main():
    parser = argparse.ArgumentParser(description="Count unique taxa by kingdom")
    parser.add_argument("input_file", help="Path to biosample_result.txt")
    args = parser.parse_args()

    organisms = parse_biosample_file(args.input_file)

    seen_taxa = {}
    kingdom_counts = defaultdict(set)
    unknown_organisms = []

    for i, organism in enumerate(set(organisms), 1):  # deduplicate by name before querying
        tax_id, kingdom = get_tax_id_and_kingdom(organism)
        if tax_id and kingdom in {"bacteria", "fungi", "virus", "other", "Unknown"}:
            if tax_id not in seen_taxa:
                seen_taxa[tax_id] = kingdom
                kingdom_counts[kingdom].add(tax_id)
        elif kingdom == "Unknown":
            unknown_organisms.append(organism)

    if unknown_organisms:
        print("\nOrganisms that could not be classified:")
    for org in unknown_organisms:
        print(f"- {org}")

        print(f"[{i}/{len(set(organisms))}] Processed: {organism} → {kingdom} (TaxID: {tax_id})")

    print("\nUnique taxa counts by kingdom:")
    print(f"fungi\t{len(kingdom_counts['fungi'])}")
    print(f"bacteria\t{len(kingdom_counts['bacteria'])}")
    print(f"virus\t{len(kingdom_counts['virus'])}")
    print(f"other\t{len(kingdom_counts['other'])}")
    print(f"unknown\t{len(kingdom_counts['Unknown'])}")

if __name__ == "__main__":
    main()
