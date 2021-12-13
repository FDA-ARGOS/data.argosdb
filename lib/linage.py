#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
from Bio import Entrez

bioproject = {}

with open('test.csv', 'r') as file:
    data = csv.reader(file)
    header = next(data)
    for row in data: 
        bioproject[row[0]] = row[1:]
print(bioproject)
test_id = bioproject['GCA_013267415.1'][5]
print(test_id)
# print(",".join(assembly_ids))
Entrez.email = "hadley_king@gwu.edu"
# print(Entrez.epost("assembly", id=",".join(assembly_ids)).read())
handle = Entrez.esearch(db="taxonomy", term=test_id)
record = Entrez.read(handle)
handle.close()
print(record['IdList'])
handle2 = Entrez.efetch(db='taxonomy', id=record['IdList'][0])
record2 = Entrez.read(handle2)
print(record2[0]['Lineage'])