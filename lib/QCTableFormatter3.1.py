# This script simply pulls data from HIVE as it is dumped and collates it into the order required by the ARGOSDB project.
# Built as a python script with functions to make it easier to manipulate in the future when changes are made to the order (i.e. just rearrange line 95).
# May be a good idea to version this script along with the data dictionary if/when that happens to keep track more easily.

import re
# Separator for when multiple values are nested in one cell:
delim = "|"
# Separator for output of the table:
sep = ';'

def codonQCReformatter(inFile):
	formattedNames = ''
	formattedNotCodList = ''
	formattedProtCodList = ''
	nameList = []
	notCodList = []
	protCodList = []

	for line in inFile:
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
			values.append(int(line.split(',')[1]))
		except:
			pass
	prctCplx = round ((values[0]/(values[0]+values[1]))*100, 1)
	prctNtCplx = round ((values[1]/(values[0]+values[1]))*100, 1)
	return prctCplx, prctNtCplx

def letterQuality(inFile):
	for line in inFile:
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

def inputGrabber(inFile):
	HIVEIDs = ''
	for line in inFile:
		HIVEIDs = HIVEIDs + str(line.rstrip())
		return str(HIVEIDs)

with open ("TEMPIDFile.txt") as IDList, open ("TEMPcodonQCTable.csv") as CQC, open ("TEMPcountNsPercentageTable.csv") as NCT, open ("TEMPComplexityTable.csv") as Complexity, open ("TEMPLetterCountQuality.csv") as LCQ:
	IDs = inputGrabber(IDList)
	print (type(IDs))
	codonName, codonCoding, codonNotCoding = codonQCReformatter(CQC)
	print (type(IDs))
	NPercent, NCount = NCountReformatter(NCT)
	print (type(IDs))
	PercentComplex, PercentNotComplex = ComplexityReformatter(Complexity)
	print (type(IDs))
	count_a, avg_quality_a, count_c, avg_quality_c, count_g, avg_quality_g, count_t, avg_quality_t = letterQuality(LCQ)
	print (type(IDs))
	print (IDs + sep + codonName + sep + codonCoding + sep + codonNotCoding + sep + NCount + sep + NPercent + sep + str(PercentComplex) + sep + str(PercentNotComplex) + sep + avg_quality_a + sep + avg_quality_t + sep + avg_quality_g + sep + avg_quality_c + sep + count_a + sep + count_c + sep + count_g + sep + count_t)


