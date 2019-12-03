import sys
import os
import urllib
import urllib2
import gzip
from math import log

GSEID = sys.argv[1]
GPLID = sys.argv[2]
geneListFile = 'data/Genes_For_EMT_Score.txt'
normalizerProbeListFile = 'data/Normalizer_Probes.txt'

GSENum = GSEID[GSEID.index('E') + 1:]
GPLNum = GPLID[GPLID.index('L') + 1:]

urlFill = ''
if len(GSENum) < 4:
	urlFill += 'nnn'
else:
	for i in range(len(GSENum) - 3):
		urlFill += GSENum[i]
	urlFill += 'nnn'
GSEurlFill = urlFill

GSEURL = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE' + urlFill + '/' + GSEID + '/matrix/' + GSEID + '_series_matrix.txt.gz'

urlFill = ''
if len(GPLNum) < 4:
	urlFill += 'nnn'
else:
	for i in range(len(GPLNum) - 3):
		urlFill += GPLNum[i]
	urlFill += 'nnn'

GPLURL = 'ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL' + urlFill + '/' + GPLID + '/annot/' + GPLID + '.annot.gz'

GPLFileName = 'annot.txt.gz'
try:
	numTries = 0
	while True:
		try:
			numTries += 1
			urllib.urlretrieve(GPLURL, GPLFileName)
			break
		except Exception, e:
			if numTries > 10:
				raise Exception('I googled your symptoms and it says you may have network connectivity issues - Andy Dwyer.')
except Exception, e:
	GPLURLalt = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + GPLID
	numTries = 0
	while True:
		try:
			numTries += 1
			GEOPage = urllib2.urlopen(GPLURLalt)
			break
		except Exception, e:
			if numTries > 10:
				raise Exception('I googled your symptoms and it says you may have network connectivity issues - Andy Dwyer.')
				sys.exit(0)
	GEOInfo = GEOPage.read()
	GEOPage.close()
	GEOLines = GEOInfo.split('\n')
	newGPLURLComp1 = ''
	newGPLURLComp2 = ''
	for line in GEOLines:
		if 'Download full table...' in line or 'View full table...' in line:
			index = line.index('&amp;id=') + len('&amp;id=')
			index2 = line.index('&amp;db=')
			indexDB = index2 + len('&amp;db=')
			indexDB2 = 0
			if 'View full table...' in line:
				indexDB2 = line.index('\', \'_blank\'')
			else:
				indexDB2 = line.index('\', \'_self\'')
			newGPLURLComp1 = line[index:index2]
			newGPLURLComp2 = line[indexDB:indexDB2]
			break
	newGPLURL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?mode=raw&is_datatable=true&acc=' + GPLID + '&id=' + newGPLURLComp1 + '&db=' + newGPLURLComp2
	GPLFileName2 = 'tempGPL.log'
	GPLFileName3 = 'tempGPL2.log'
	urllib.urlretrieve(newGPLURL, GPLFileName2)
	with open(GPLFileName2, 'r') as f, open(GPLFileName3, 'w') as g:
		for line in f:
			if line[0] == '#':
				continue
			if 'ID' in line:
				l = line.strip().split('\t')
				for val in l:
					if 'SYMBOL' in val.upper():
						g.write('Gene symbol\t')
					else:
						g.write(val.strip() + '\t')
				g.write('\n')
			else:
				g.write(line.strip() + '\n')
	fIn = open(GPLFileName3, 'r')
	fOut = gzip.open(GPLFileName, 'w')
	fOut.writelines(fIn)
	fIn.close()
	fOut.close()

	os.remove(GPLFileName2)
	os.remove(GPLFileName3)

geneList = []
with open(geneListFile, 'r') as f:
	for line in f:
		l = line.strip()
		if len(l) > 0:
			geneList.append(l)

GeneIDMap = {}
for gene in geneList:
	GeneIDMap[gene] = []

flag = 0
probeIndex = 0
nameIndex = 0
with gzip.open(GPLFileName, 'r') as f:
	for line in f:
		l = line.strip().split('\t')
		if len(l) > 3 and flag == 1:
			probeID = l[probeIndex].strip()
			geneName = l[nameIndex].strip().split('//')
			for kindex in range(len(geneName)):
				if '/' in geneName[kindex]:
					geneName[kindex] = geneName[kindex][1:]
			for gene in geneName:
				if gene.strip().upper() in geneList:
					GeneIDMap[gene.strip().upper()].append(probeID)
		if 'ID' in l and ('Gene symbol' in l or 'gene_assignment' in l or 'ORF' in l):
			probeIndex = l.index('ID')
			nameIndex = -1
			if 'Gene symbol' not in l and 'ORF' not in l:
				nameIndex = l.index('gene_assignment')
			elif 'gene_assignment' not in l and 'ORF' not in l:
				nameIndex = l.index('Gene symbol')
			else:
				nameIndex = l.index('ORF')
			flag = 1

GSEFileName = 'raw_data.txt.gz'
try:
	numTries = 0
	while True:
		try:
			numTries += 1
			urllib.urlretrieve(GSEURL, GSEFileName)
			break
		except Exception, e:
			if numTries > 10:
				raise Exception('I googled your symptoms and it says you may have network connectivity issues - Andy Dwyer.')
except Exception, e:
	GSEURLIDalt = GSEID + '-' + GPLID
	GSEURLalt = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE' + GSEurlFill + '/' + GSEID + '/matrix/' + GSEURLIDalt + '_series_matrix.txt.gz'
	numTries = 0
	while True:
		try:
			numTries += 1
			urllib.urlretrieve(GSEURLalt, GSEFileName)
			break
		except Exception, e:
			if numTries > 10:
				raise Exception('I googled your symptoms and it says you may have network connectivity issues - Andy Dwyer.')
				sys.exit(0)

takeLogFlag = 0

patientNum = 0
flag = 0

with gzip.open(GSEFileName, 'r') as f:
	for line in f:
		if flag == 1:
			l = line.strip().split('\t')
			if len(l) != patientNum:
				continue
			flag2 = 0
			for val in l:
				if len(val) < 1 or 'null' in val.strip().lower():
					flag2 = 1
			if flag2 == 1:
				continue
			probe = l[0].strip()
			for val in l[1:]:
				expressionValue = float(val.strip())
				if expressionValue > 100.0:
					takeLogFlag = 1
		if 'ID_REF' in line:
			l = line.strip().split('\t')
			patientNum = len(l)
			flag = 1

GSEOutFile = open(GSEID + 'EMTScoreData.log', 'w')

GSEData = {}
probeCounts = {}
for gene in geneList:
	probeCounts[gene] = 0

patientNum = 0
flag = 0

with gzip.open(GSEFileName, 'r') as f:
	for line in f:
		if flag == 1:
			l = line.strip().split('\t')
			if len(l) != patientNum:
				continue
			flag2 = 0
			for val in l:
				if len(val) < 1 or 'null' in val.strip().lower():
					flag2 = 1
			if flag2 == 1:
				continue
			probe = l[0].strip()
			if '"' in probe:
				probe = probe[1:-1]
			for gene in geneList:
				if probe in GeneIDMap[gene]:
					probeCounts[gene] += 1
					probeData = []
					for val in l[1:]:
						if 'null' in val:
							probeData.append(0.0)
						else:
							expressionValue = float(val.strip())
							if takeLogFlag == 1:
								expressionValue = log(expressionValue, 2)
							probeData.append(expressionValue)
					if gene in GSEData and len(GSEData[gene]) > 0:
						for i in range(len(probeData)):
							GSEData[gene][i] += probeData[i]
					else:
						GSEData[gene] = []
						for val in probeData:
							GSEData[gene].append(val)
		if 'ID_REF' in line:
			l = line.strip().split('\t')
			for val in l:
				if '"' in l:
					GSEOutFile.write(val.strip()[1:-1] + ' ')
				else:
					GSEOutFile.write(val.strip() + ' ')
			GSEOutFile.write('\n')
			patientNum = len(l)
			flag = 1

for gene in GSEData:
	for i in range(len(GSEData[gene])):
		GSEData[gene][i] /= probeCounts[gene]

for gene in geneList[0:5]:
	GSEOutFile.write(gene + ' ')
	for val in GSEData[gene]:
		GSEOutFile.write(str(val) + ' ')
	GSEOutFile.write('\n')

normalizerProbes = []
with open(normalizerProbeListFile, 'r') as f:
	for line in f:
		l = line.strip()
		normalizerProbes.append(l)

availProbes = []
with open('data/Available_NCI60_Probes.txt', 'r') as f:
	for line in f:
		l = line.strip()
		availProbes.append(l)

GSENormalizerData = {}
patientNum = 0
flag = 0

with gzip.open(GSEFileName, 'r') as f:
	for line in f:
		if flag == 1:
			l = line.strip().split('\t')
			if len(l) != patientNum:
				continue
			flag2 = 0
			for val in l:
				if len(val) < 1 or 'null' in val.strip().lower():
					flag2 = 1
			if flag2 == 1:
				continue
			probe = l[0].strip()
			if '"' in probe:
				probe = probe[1:-1]
			if probe in normalizerProbes and probe in availProbes:
				GSENormalizerData[probe] = []
				for val in l[1:]:
					if 'null' in val:
						GSENormalizerData[probe].append(0.0)
					else:
						expressionValue = float(val.strip())
						if takeLogFlag == 1:
							expressionValue = log(expressionValue, 2)
						GSENormalizerData[probe].append(expressionValue)
		if 'ID_REF' in line:
			l = line.strip().split('\t')
			patientNum = len(l)
			flag = 1

usedNormalizers = []
usedNormalizersProbeIndices = []
for probe in normalizerProbes:
	if probe not in GSENormalizerData:
		continue
	GSEOutFile.write(probe + ' ')
	for val in GSENormalizerData[probe]:
		GSEOutFile.write(str(val) + ' ')
	GSEOutFile.write('\n')
	usedNormalizers.append(probe)
	usedNormalizersProbeIndices.append(availProbes.index(probe))

with open(GSEID + 'IndicesUsedNormalizers.log', 'w') as f:
	for index in usedNormalizersProbeIndices:
		f.write(str(index + 1) + '\n')

os.remove(GSEFileName)
os.remove(GPLFileName)
