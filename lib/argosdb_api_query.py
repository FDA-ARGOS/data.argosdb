#!/usr/bin/env python3
"""Argos DB API Query

Entry point for ARGOS DB
"""

__version__ = "0.1.0"
__status__ = "BETA"

import os
import requests
api_url = "https://api.argosdb.org/records/search"
results = []
assemblies = []
data = [{
  "bcoid": "ARGOS_000012",
  "offset": 1,
  "limit": 10000
},
{
  "bcoid": "ARGOS_000022",
  "offset": 1,
  "limit": 10000
}]

os.system('mkdir home/assembly')

for item in data:
    response = requests.post(api_url, json=item)
    results.append(response.json())
    for record in response.json()['recordlist']:
        if record['genome_assembly_id'] not in assemblies:
            assemblies.append(record['genome_assembly_id'])
            os.system(f"efetch -db assembly -id {record['genome_assembly_id']} -format docsum > home/{record['genome_assembly_id']}.xml")
    print(response.status_code)
print(assemblies)
