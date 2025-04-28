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

# def parse_biosample_file(filepath):
#     with open(filepath, 'r') as f:
#         content = f.read()

#     entries = re.split(r'\n(?=\d+\. Pathogen:)', content)
#     #entries = re.split(r'(?=\d+\. Pathogen:)', content)
#     organisms = []

#     # for entry in entries:
#     #     #print(entries)
#     #     organism_match = re.search(r'Organism:\s*(.+)', entry)
#     #     if organism_match:
#     #         organisms.append(organism_match.group(1).strip())
#     #         #print(f" Organism: {organism_match.group(1).strip()}")
#     #     else:
#     #         print(f"⚠️ Warning: No organism found in entry \n {idx}\n")

#     for idx, entry in enumerate(entries, 1):
#         organism_match = re.search(r'Organism:\s*(.+)', entry)
#         if organism_match:
#             organisms.append(organism_match.group(1).strip())
#         else:
#             print(f"⚠️ Warning: No organism found in entry {idx}")

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

    print(f"✅ Parsed {len(organisms)} organisms from file.")
    return organisms


def get_kingdom(organism):
    api_key = os.getenv("NCBI_API_KEY") or "YOURKEY HERE"                        #<--- add your api key here

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
        print('organism failed stat code:', organism)
        return "Unknown"

    try:
        tax_id = search_resp.json()["esearchresult"]["idlist"][0]
    except Exception:
        print("no idlist found---", organism, "\n")
        return "Unknown"

    # Fetch lineage
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
        print("bad status code------------", organism)
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
            print('--- for Other org:   ', lineage)
            return "Other"
    except Exception:
        return "Unknown"

def main():
    parser = argparse.ArgumentParser(description="Count BioSamples by Kingdom")
    parser.add_argument("input_file", help="Path to biosample_result.txt")
    args = parser.parse_args()

    organisms = parse_biosample_file(args.input_file)

    kingdom_counts = Counter()
    unknown_organisms = []

    for i, organism in enumerate(organisms, 1):
        kingdom = get_kingdom(organism)
        kingdom_counts[kingdom] += 1
        if kingdom == "Unknown":
            unknown_organisms.append(organism)
        #print(f"[{i}/{len(organisms)}] {organism} => {kingdom}")

    if unknown_organisms:
        print("\nOrganisms that could not be classified:")
    for org in unknown_organisms:
        print(f"- {org}")


    print("\nBiosample counts by kingdom:")
    print(f"fungi\t{kingdom_counts.get('fungi', 0)}")
    print(f"bacteria\t{kingdom_counts.get('bacteria', 0)}")
    print(f"virus\t{kingdom_counts.get('virus', 0)}")
    print(f"other\t{kingdom_counts.get('Other', 0)}")
    print(f"unknown\t{kingdom_counts.get('Unknown', 0)}")
    print(f"total: \t{sum(kingdom_counts.values())}")
    print(f"expected total organisms: {len(organisms)}")

    # if sum(kingdom_counts.values()) != len(organisms):
    #     print("\n⚠️ Warning: Counts do not add up! Some entries may have been missed.")

if __name__ == "__main__":
    main()
