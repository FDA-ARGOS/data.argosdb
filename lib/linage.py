#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv, time
from Bio import Entrez
csv_file = 'data_files/PRJNA231221_AssemblyDetails.csv'
Entrez.email = "hadley_king@gwu.edu"
bioproject = {}

with open(csv_file, 'r') as file:
    data = csv.reader(file)
    header = next(data)
    for row in data: 
        bioproject[row[0]] = row[1:]

header.append('Taxonomy ID')
header.append('Full Lineage')
with open('csv_file.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)
    for key in bioproject:
        handle = Entrez.esearch(db="taxonomy", term=bioproject[key][5])
        record = Entrez.read(handle)
        handle.close()
        print(record)
        try: 
            bioproject[key].append(record['IdList'][0])
            print(record['IdList'])
            handle2 = Entrez.efetch(db='taxonomy', id=record['IdList'][0])
            record2 = Entrez.read(handle2)
            print(record2[0]['Lineage'])
            bioproject[key].append(record2[0]['Lineage'])
        except: 
            bioproject[key].append('')
            bioproject[key].append('')
        time.sleep(1)
        row = [key]
        for item in bioproject[key]:
            row.append(item)
        writer.writerow(row)
# print(",".join(assembly_ids))
# print(Entrez.epost("assembly", id=",".join(assembly_ids)).read())

