# christie woodside
# april 28, 2025

'''counts the assemblies for each organism for each kingdom for the NCBI bioproject. takes in 
the assembly file which is the list of the assemblies from the ncbi bioproject itself.  This can be found when you click on the assembly hyperlink
from the bioproject page > download table from the assembly dataset page'''



import pandas as pd
import requests
import os
import time
import xml.etree.ElementTree as ET
from collections import Counter

def get_kingdom_from_taxid(tax_id):
    api_key = os.getenv("NCBI_API_KEY") or "YOUR API KEY HERE"    <------ add your NCBI api key
    
    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "taxonomy",
        "id": tax_id,
        "retmode": "xml",
        "api_key": api_key
    }
    
    response = requests.get(efetch_url, params=params)
    if response.status_code == 429:
        time.sleep(1)
        response = requests.get(efetch_url, params=params)
    if response.status_code != 200:
        print(f"⚠️ Failed for tax_id {tax_id}")
        return "Unknown"
    
    try:
        root = ET.fromstring(response.content)
        lineage = root.findtext("Taxon/Lineage")
        if lineage:
            if "Viruses" in lineage:
                return "virus"
            elif "Fungi" in lineage:
                return "fungi"
            elif "Bacteria" in lineage:
                return "bacteria"
            else:
                print(f"--- Other lineage for tax_id {tax_id}: {lineage}")
                return "Other"
        else:
            return "Unknown"
    except Exception:
        return "Unknown"

def main():
    input_file = "/Users/user/Desktop/assembly_result.tsv"                # <<< UPDATE this to your filepath
    df = pd.read_csv(input_file, sep='\t')
    df.columns = df.columns.str.strip()

    # Drop duplicates based on Assembly Name
    df_unique = df.drop_duplicates(subset=['Assembly Name'])

    kingdom_counts = Counter()
    unknown_taxids = []
    taxid_cache = {}  # optional small cache to avoid redundant queries

    #print(df_unique.columns.tolist())
    # exit()


    for idx, (i, row) in enumerate(df_unique.iterrows(), 1):  # <= change to iterrows()
        tax_id = row['Organism Taxonomic ID']
        assembly_name = row['Assembly Name']

        if pd.isna(tax_id):
            print(f"⚠️ Missing tax_id for {assembly_name}")
            kingdom = "Unknown"
        else:
            try:
                tax_id_str = str(int(tax_id))  # make sure it’s a clean string
            except ValueError:
                print(f"⚠️ Skipping invalid tax_id: {tax_id}")
                kingdom = "Unknown"
            else:
                # Check cache
                if tax_id_str in taxid_cache:
                    kingdom = taxid_cache[tax_id_str]
                else:
                    kingdom = get_kingdom_from_taxid(tax_id_str)
                    taxid_cache[tax_id_str] = kingdom

        kingdom_counts[kingdom] += 1
        if kingdom == "Unknown":
            unknown_taxids.append(tax_id)

        # if idx % 50 == 0 or idx == len(df_unique):
        #     print(f"[{idx}/{len(df_unique)}] {assembly_name} => {kingdom}")


    print("\nAssembly counts by kingdom:")
    print(f"fungi\t{kingdom_counts.get('fungi', 0)}")
    print(f"bacteria\t{kingdom_counts.get('bacteria', 0)}")
    print(f"virus\t{kingdom_counts.get('virus', 0)}")
    print(f"other\t{kingdom_counts.get('Other', 0)}")
    print(f"unknown\t{kingdom_counts.get('Unknown', 0)}")
    print(f"total: \t{sum(kingdom_counts.values())}")
    print(f"expected total assemblies: {len(df_unique)}")

    if unknown_taxids:
        print("\nTax IDs that could not be classified:")
        for tax in unknown_taxids:
            print(f"- {tax}")

if __name__ == "__main__":
    main()
