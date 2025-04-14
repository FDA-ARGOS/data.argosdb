import re
import requests
import csv
import os
import argparse
import time
import xml.etree.ElementTree as ET

def parse_biosample_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    entries = re.split(r'\n(?=\d+\. Pathogen:)', content)
    samples = []

    for entry in entries:
        biosample_match = re.search(r'BioSample:\s*(SAMN\d+)', entry)
        strain_match = re.search(r'strain:\s*(\S+)', entry)
        organism_match = re.search(r'Organism:\s*(.+)', entry)

        if biosample_match:
            biosample_id = biosample_match.group(1)
            strain = strain_match.group(1) if strain_match else ""
            organism = organism_match.group(1).strip() if organism_match else ""
            samples.append({
                "biosample_id": biosample_id,
                "strain": strain,
                "organism": organism
            })

    return samples

def is_virus(biosample_id):
    api_key = os.getenv("NCBI_API_KEY") or "bfbde99c962d228023e8d62a078bdb12d108"

    # Step 1: Use esearch to get the UID from the BioSample accession
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    esearch_params = {
        "db": "biosample",
        "term": biosample_id,
        "retmode": "json",
        "api_key": api_key
    }
    esearch_resp = requests.get(esearch_url, params=esearch_params)
    if esearch_resp.status_code == 429:
        time.sleep(1)
        esearch_resp = requests.get(esearch_url, params=esearch_params)
    if esearch_resp.status_code != 200:
        print(f"Warning: Failed to esearch {biosample_id} (status {esearch_resp.status_code})")
        return False

    try:
        uid_list = esearch_resp.json().get("esearchresult", {}).get("idlist", [])
        if not uid_list:
            print(f"Warning: No UID found for {biosample_id}")
            return False
        uid = uid_list[0]
    except Exception as e:
        print(f"Error extracting UID for {biosample_id}: {e}")
        return False

    # Step 2: Get the BioSample summary using the UID
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    summary_params = {
        "db": "biosample",
        "id": uid,
        "retmode": "json",
        "api_key": api_key
    }
    summary_response = requests.get(summary_url, params=summary_params)
    if summary_response.status_code == 429:
        time.sleep(1)
        summary_response = requests.get(summary_url, params=summary_params)
    if summary_response.status_code != 200:
        print(f"Warning: Failed to retrieve summary for UID {uid} (status {summary_response.status_code})")
        return False

    try:
        data = summary_response.json()
        doc = data['result'][uid]
        organism = doc.get('organism', '')
    except Exception as e:
        print(f"Warning: No summary data for {biosample_id} UID {uid}: {e}")
        return False

    if not organism:
        print(f"Warning: No organism found for {biosample_id}")
        return False

    # Step 3: Get taxonomy ID from organism name
    tax_search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    tax_params = {
        "db": "taxonomy",
        "term": organism,
        "retmode": "json",
        "api_key": api_key
    }
    tax_search_resp = requests.get(tax_search_url, params=tax_params)
    if tax_search_resp.status_code == 429:
        time.sleep(1)
        tax_search_resp = requests.get(tax_search_url, params=tax_params)
    if tax_search_resp.status_code != 200:
        print(f"Warning: Failed to search taxonomy for {organism} (status {tax_search_resp.status_code})")
        return False

    try:
        tax_id_list = tax_search_resp.json().get("esearchresult", {}).get("idlist", [])
        if not tax_id_list:
            print(f"Warning: No taxonomy ID found for {organism}")
            return False
        tax_id = tax_id_list[0]
    except Exception as e:
        print(f"Error parsing taxonomy ID for {organism}: {e}")
        return False

    # Step 4: Get full lineage using efetch and check for 'Viruses'
    tax_fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    tax_fetch_params = {
        "db": "taxonomy",
        "id": tax_id,
        "retmode": "xml",
        "api_key": api_key
    }
    tax_fetch_resp = requests.get(tax_fetch_url, params=tax_fetch_params)
    if tax_fetch_resp.status_code == 429:
        time.sleep(1)
        tax_fetch_resp = requests.get(tax_fetch_url, params=tax_fetch_params)
    if tax_fetch_resp.status_code != 200:
        print(f"Warning: Failed to fetch taxonomy XML for tax ID {tax_id} (status {tax_fetch_resp.status_code})")
        return False

    try:
        root = ET.fromstring(tax_fetch_resp.content)
        lineage = root.findtext("Taxon/Lineage")
        return lineage and "Viruses" in lineage
    except Exception as e:
        print(f"Error parsing XML lineage for tax ID {tax_id}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Filter virus samples from BioSample file")
    parser.add_argument("input_file", help="Path to input biosample_result.txt file")
    parser.add_argument("output_file", help="Path to output TSV file")
    args = parser.parse_args()

    samples = parse_biosample_file(args.input_file)
    virus_samples = []

    for sample in samples:
        biosample_id = sample['biosample_id']
        if is_virus(biosample_id):
            virus_samples.append({
                "Organism": sample["organism"],
                "BioSample ID": biosample_id,
                "Strain": sample["strain"]
            })

    with open(args.output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Organism", "BioSample ID", "Strain"], delimiter='\t')
        writer.writeheader()
        writer.writerows(virus_samples)

    print(f"Done. Found {len(virus_samples)} virus samples. Output written to {args.output_file}")

if __name__ == "__main__":
    main()