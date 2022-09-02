#%%
#!/usr/bin/env python3
"""assemblyQC Data Grabber

Designed to grab and modify data the number of unidentified nucleotides (num_N) from HIVE NGS Multi-QC processes
"""
import requests
import getpass
import json
import re
import csv
import os

# Separator for when multiple values are nested in one cell:
delim = "|"
# Separator for output of the table:
sep = ';'

#_______________________________________________________________# Functions for formatting the data.

def codon_qc_reformatter(input_file):
    """Codon QC Reformatter
    
    IDK???

    Parameters
    ----------
    input_file: str

    """

    formattedNames = ''
    formattedNotCodList = ''
    formattedProtCodList = ''
    nameList = []
    notCodList = []
    protCodList = []

    for line in input_file:
        name = re.search("\".*\"", line)
        coding = re.split("\".*\"", line.strip())
        notCod = int(coding[1].split(",")[1])
        protCod = int(coding[1].split(",")[2])
        totCod = notCod+protCod
        nameList.append(name.group(0)[1:-1])
        notCodList.append(round(notCod/totCod, 2))
        protCodList.append(round(protCod/totCod, 2))

    for i in nameList:
        formattedNames = formattedNames + str(i) + delim
    for j in notCodList:
        formattedNotCodList = formattedNotCodList + str(j) + delim
    for k in protCodList:
        formattedProtCodList = formattedProtCodList + str(k) + delim
    return formattedNames[:-1], formattedProtCodList[:-1], formattedNotCodList[:-1]

def NCountReformatter(Ns):
    percentage = []
    count = []
    NPercent = ''
    NCount = ''
    for line in Ns:
        percentage.append(line.split(',')[0])
        count.append(line.split(',')[1])
    for i in percentage:
        NPercent = NPercent + str(i) + delim
    for j in count:
        NCount = NCount + str(j).rstrip() + delim
    return NPercent[:-1], NCount[:-1]

def ComplexityReformatter(Cmplx):
    values = []
    for line in Cmplx:
        try:
#                       int(line.split(',')[1])
            values.append(int(line.split(',')[1]))
        except:
            pass
    prctCplx = round ((values[0]/(values[0]+values[1]))*100, 1)
    prctNtCplx = round ((values[1]/(values[0]+values[1]))*100, 1)
    return prctCplx, prctNtCplx

def letterQuality(input_file):
    for line in input_file:
        if line.split(',')[0] == 'A':
            ct_a = line.split(',')[1]
            qual_a = line.split(',')[2]
        elif line.split(',')[0] == 'C':
            ct_c = line.split(',')[1]
            qual_c = line.split(',')[2]
        if line.split(',')[0] == 'G':
            ct_g = line.split(',')[1]
            qual_g = line.split(',')[2]
        if line.split(',')[0] == 'T':
            ct_t = line.split(',')[1]
            qual_t = line.split(',')[2]
        else:
            pass
    return ct_a, qual_a.rstrip(), ct_c, qual_c.rstrip(), ct_g, qual_g.rstrip(), ct_t, qual_t.rstrip()

def inputGrabber(input_file):
    HIVEIDs = ''
    for line in input_file:
        HIVEIDs = HIVEIDs + str(line.rstrip())
        return str(HIVEIDs)

#_______________________________________________________________# Pick HIVE instance.

u1 = "https://hive.biochemistry.gwu.edu/dna.cgi?"
u2 = "https://hive2.biochemistry.gwu.edu/dna.cgi?"
ua = "https://hive.aws.biochemistry.gwu.edu/dna.cgi?"

userHIVE = None

def getHIVE():
    global userHIVE
    userHIVE = input("\nAre you using HIVE1, HIVE2, or HIVE AWS (type 1, 2, or a and hit enter)? ")

getHIVE()

if userHIVE == '1':
    basename = u1
    hiveVal = 'HIVE 1'
elif userHIVE == '2':
    basename = u2
    hiveVal = 'HIVE 2'
elif userHIVE == 'a':
    basename = ua
    hiveVal = 'HIVE AWS'
else:
    print ('wrong input detected, please try again.')
    getHIVE()

print ("\nUsing " + basename[0:-8] + "... ")

#_______________________________________________________________# Login.

user = input("\nEnter your username on " + hiveVal + " (email address): ")
password  = getpass.getpass()
loginParams = {'api': '0', 'cmdr': 'login', 'login': user, 'pswd': password}
response = requests.get(basename, params=loginParams)
cookies = response.cookies

"""
Maybe consider putting something here that will only advance if login was successful. E.g. if response == 200, else throw login error.
"""

#_______________________________________________________________# Get Object ID.

userVal = input ("Enter the Object ID for the Multiple NGS analysis that you'd like to pull data from: ")
print ("\nAttempting to retrieve data from object " + userVal)
print ("\n")

params1 = {'cmdr': '-qpRawSubmit', 'check': '1', 'svc': 'tblqryx4', 'oper': 'list', 'cnt': '-1', 'raw': '1', 'cols': '0-20', 'objs': str(userVal), 'tbl': 'mqc.fileTable.csv'}
params2 = {'cmdr': '-qpRawSubmit', 'check': '1', 'svc': 'tblqryx4', 'oper': 'list', 'cnt': '-1', 'raw': '1', 'cols': '0-20', 'objs': str(userVal), 'tbl': 'mqc.positionTable.csv'}

response1 = requests.get(basename, params1, cookies=cookies)
response2 = requests.get(basename, params2, cookies=cookies)

reqGrab1 = response1.text.splitlines()
reqGrab2 = response2.text.splitlines()

reqID1 = reqGrab1[1].split(",")[4]
reqID2 = reqGrab2[1].split(",")[4]

params3 = {'cmdr': '-qpData', 'req': reqID1, 'raw': '1', 'grp': '1', 'dname': '_.csv'}
params4 = {'cmdr': '-qpData', 'req': reqID2, 'raw': '1', 'grp': '1', 'dname': '_.csv'}
response3 = requests.get(basename, params3, cookies=cookies)
response4 = requests.get(basename, params4, cookies=cookies)

with open ("mqc.fileTable.csv", "w+") as outFile1:
    print (response3.text, file=outFile1)

with open ("mqc.positionTable.csv", "w+") as outFile2:
    print (response4.text, file=outFile2)

#_______________________________________________________________# Grab input read IDs.

inputParams = {'cmdr': 'propget', 'mode': 'json', 'ids': userVal, 'files': '*'}
responseInput = requests.get(basename, inputParams, cookies=cookies)

lst =[]
meta = json.loads(responseInput.text)
data = meta["query_objId"]
print(meta)
for key, value in data.items():
        lst.append(value)

#_______________________________________________________________# Get SRR/ERR archive accession for each file.

try:
    print ("Getting data from input read files...\n")
except:
    print ("Error: cannot find read files for the computation ID you entered. Did you login correctly and enter the correct \"Multiple ngs QC\" Object ID corresponding to the instance of HIVE you chose?\n")
    exit()

strIDs = ','.join(str(v) for v in lst)
SRRParams = {'cmdr': 'propget', 'prop': 'id,_brief', 'mode': 'csv', 'ids': strIDs, 'raw': '1'}
responseSRR = requests.get(basename, SRRParams, cookies=cookies)
decoded_content = responseSRR.content.decode('utf-8')
cr = csv.reader(decoded_content.splitlines(), delimiter=',')
with open('TEMPIDFile.csv', 'w+', newline='') as IDFile:
    my_list = list(cr)
    for row in my_list:
        print(row[3].split()[0], file=IDFile)

#_______________________________________________________________# Pull data files for each input read.

"""
These API calls are basically repeated 4 times. They should just be pulled into a function.
"""

with open ("outputFile.csv", "w+") as finalFile:
    print ("file_name;codon_table;%protein_coding;%not_coding;%_reads;density_Ns_per_read;%complex;%not_complex;avg_quality_a;avg_quality_t;avg_quality_g;avg_quality_c;count_a;count_c;count_g;count_t", file=finalFile)
    with open ('TEMPIDFile.csv') as IDFile:
        next(iter(IDFile))
        for a,b in zip(lst,IDFile):
            print ("HIVE Object ID: " + str(a))
            print ("Input file: " + str(b).rstrip())
            with requests.Session() as s:
                with open ("TEMPcodonQCTable.csv", "w+") as outFileCodon:
                    codonParams = {'cmdr': 'objFile', 'filename': '.qc2.codonQCTable.csv', 'ids': a,}
                    codonInput = s.get(basename, params=codonParams, cookies=cookies)
                    codonDecode = codonInput.content.decode('utf-8')
                    codD = csv.reader(codonDecode.splitlines(), delimiter=',')
                    codonList = list(codD)
                    for row in codonList:
                        print ("\"" + row[0]+ "\"" + ','+row[1]+','+row[2], file = outFileCodon)

                with open ("TEMPLetterCountQuality.csv", "w+") as outFileLetter:
                    letterParams = {'cmdr': 'objFile', 'filename': '.qc2.sumLetterTable.csv', 'ids': a,}
                    letterInput = s.get(basename, params=letterParams, cookies=cookies)
                    letterDecode = letterInput.content.decode('utf-8')
                    letD = csv.reader(letterDecode.splitlines(), delimiter=',')
                    letterList = list(letD)
                    for row in letterList:
                        print (row[0]+','+row[1]+','+row[2], file = outFileLetter)

                with open ("TEMPcountNsPercentageTable.csv", "w+") as outFileCount:
                    countParams = {'cmdr': 'objFile', 'filename': '.qc2.countNsPercentageTable.csv', 'ids': a,}
                    countInput = requests.get(basename, params=countParams, cookies=cookies)
                    countDecode = countInput.content.decode('utf-8')
                    couD = csv.reader(countDecode.splitlines(), delimiter=',')
                    countList = list(couD)
                    for row in countList:
                        print (row[0]+','+row[1], file = outFileCount)

                with open ("TEMPComplexityTable.csv", "w+") as outFileComplex:
                    complexParams = {'cmdr': 'objFile', 'filename': '.qc2.ComplexityTable.csv', 'ids': a,}
                    complexInput = requests.get(basename, params=complexParams, cookies=cookies)
                    complexDecode = complexInput.content.decode('utf-8')
                    comD = csv.reader(complexDecode.splitlines(), delimiter=',')
                    complexList = list(comD)
                    for row in complexList:
                        print (row[0]+','+row[1], file = outFileComplex)

#_______________________________________________________________# Data formatter.

                with open ("TEMPcodonQCTable.csv") as CQC, open ("TEMPcountNsPercentageTable.csv") as NCT, open ("TEMPComplexityTable.csv") as Complexity, open ("TEMPLetterCountQuality.csv") as LCQ:
                    next (CQC)
                    next (NCT)
                    next (Complexity)
                    next (LCQ)
                    IDs = b.rstrip()
                    codonName, codonCoding, codonNotCoding = codon_qc_reformatter(CQC)
                    NPercent, NCount = NCountReformatter(NCT)
                    PercentComplex, PercentNotComplex = ComplexityReformatter(Complexity)
                    count_a, avg_quality_a, count_c, avg_quality_c, count_g, avg_quality_g, count_t, avg_quality_t = letterQuality(LCQ)
                    print (IDs + sep+ codonName + sep + codonCoding + sep + codonNotCoding + sep + NCount + sep + NPercent + sep + str(PercentComplex) + sep + str(PercentNotComplex) + sep + avg_quality_a + sep + avg_quality_t + sep + avg_quality_g + sep + avg_quality_c + sep + count_a + sep + count_c + sep + count_g + sep + count_t, file=finalFile)

#_______________________________________________________________# Clean up.

os.remove("TEMPcodonQCTable.csv")
os.remove("TEMPcountNsPercentageTable.csv")
os.remove("TEMPComplexityTable.csv")

with open ('TEMPLetterCountQuality.csv') as LC:
    for line in LC:
        if line[0] == 'N':
            print(line)

os.remove("TEMPLetterCountQuality.csv")

print ("\nDone.")

# %%
