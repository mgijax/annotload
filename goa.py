#!/usr/local/bin/python

'''
#
# goa.py
#
# Report:
#	TR 7904/7926
#
#	Takes the GOA file ${GOAINPUTFILE} and generates
#
#		mgi.error
#			file of GOA annotations that originated from MGI
#
#		unresolvedA.error
#			file of GOA annotations where UniProtID resolves to more than one MGI Marker
#	
#		unresolvedB.error
#			file of GOA annotations where UniProtID don't resolve to a MGI Marker
#
#		duplicates.error
#		   	file of GOA annotations that are duplicates to those in MGI
#
#		pubmed.error
#			file of GOA annotations where PubMed IDs are not in MGI
#
#		goa.mgi
#			file of GOA annotations that can be appended to MGI GO association file
#
#		goa.annot
#			file of GOA annotations that can be loaded into MGI via the Annotation loader
#
# The annotation loader format has the following columns:
#
#	A tab-delimited file in the format:
#		field 1: Accession ID of Vocabulary Term being Annotated to
#		field 2: ID of MGI Object being Annotated (ex. MGI ID)
#		field 3: J: (J:#####)
#		field 4: Evidence Code Abbreviation (max length 5)
#		field 5: Inferred From (max length 255)
#		field 6: Qualifier (max length 255)
#		field 7: Editor (max length 30)
#		field 8: Date (MM/DD/YYYY)
#		field 9: Notes (max length 255)
#
# Usage:
#       goa.py
#
# History:
#
# lec   10/02/2006
#       - created
#
'''

import sys
import os
import string
import re
import db
import mgi_utils
import reportlib
import loadlib

inFileName = os.environ['GOAINPUTFILE']
createdBy = os.environ['GOAEDITOR']

unresolvedAErrorFile = reportlib.init('unresolvedA', outputdir = os.environ['GOADIR'], printHeading = 0, fileExt = '.error')
unresolvedBErrorFile = reportlib.init('unresolvedB', outputdir = os.environ['GOADIR'], printHeading = 0, fileExt = '.error')
mgiErrorFile = reportlib.init('mgi', outputdir = os.environ['GOADIR'], printHeading = 0, fileExt = '.error')
pubmedErrorFile = reportlib.init('pubmed', outputdir = os.environ['GOADIR'], printHeading = 0, fileExt = '.error')
pubmedISSErrorFile = reportlib.init('pubmedISS', outputdir = os.environ['GOADIR'], printHeading = 0, fileExt = '.error')
dupErrorFile = reportlib.init('duplicates', outputdir = os.environ['GOADIR'], printHeading = 0, fileExt = '.error')
mgiFile = reportlib.init('goa', outputdir = os.environ['GOADIR'], printHeading = 0, fileExt = '.mgi')
annotFile = reportlib.init('goa', outputdir = os.environ['GOADIR'], printHeading = 0, fileExt = '.annot')

assoc = {}	# dictionary of GOA ID:Marker MGI ID
marker = {}	# dictionary of MGI Marker ID:Marker data
annot = {}	# list of existing Marker key, GO ID, Evidence Code, Pub Med ID annotations
annotByGOID = []
annotByRef = []
pubmed = {}	# dictionary of pubmed->J:

mgiLine = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'
annotLine = '%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t\n'

loaddate = loadlib.loaddate

#
# Mouse Markers
#

db.sql('select mgiID = a.accID, m._Marker_key, m.symbol, m.name, markerType = t.name ' + \
	'into #markers ' + \
	'from ACC_Accession a, MRK_Marker m, MRK_Types t ' + \
	'where m._Organism_key = 1 ' + \
	'and m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1 ' + \
	'and m._Marker_Type_key = t._Marker_Type_key', None)
db.sql('create index idx1 on #markers(_Marker_key)', None)

results = db.sql('select * from #markers', 'auto')
for r in results:
    marker[r['mgiID']] = r

#
# Mouse Markers annotated to...
#
# SwissProt (13)
# TrEMBL (41)
# RefSeq (27)
# ENSEMBL (60)
# VEGA (85)
#

results = db.sql('select m._Marker_key, m.mgiID, goaID = a.accID ' + \
	'from #markers m, ACC_Accession a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._LogicalDB_key in (13, 41, 27, 60, 85) ' + \
	'and a._MGIType_key = 2 ', 'auto')
for r in results:
    key = r['goaID']
    value = r['mgiID']
    if not assoc.has_key(key):
	assoc[key] = []
    assoc[key].append(value)

#
# existing GO annotations that have pub med ids
# to detect duplicate annotations
#

results = db.sql('select goID = a.accID, t._Object_key, ec.abbreviation, refID = "PMID:" + r.accID ' + \
	'from VOC_Annot t, ACC_Accession a, VOC_Evidence e, VOC_Term ec, ACC_Accession r ' + \
	'where t._AnnotType_key = 1000 ' + \
	'and t._Term_key = a._Object_key ' + \
	'and a._MGIType_key = 13 ' + \
	'and a.preferred = 1 ' + \
	'and t._Annot_key = e._Annot_key ' + \
	'and e._EvidenceTerm_key = ec._Term_key ' + \
	'and e._Refs_key = r._Object_key ' + \
	'and r._MGIType_key = 1 ' + \
	'and r._LogicalDB_key = 29', 'auto')
for r in results:

    key = r['_Object_key']

    if key not in annot:
	annot[key] = []
    annot[key].append((r['goID'], r['abbreviation'], r['refID']))

    if r['goID'] not in annotByGOID:
        annotByGOID.append(r['goID'])

    if r['refID'] not in annotByRef:
        annotByRef.append(r['refID'])

#
# existing IEA GO annotations
# J:72247 interpro
# J:60000 swissprot
# J:72245
#

results = db.sql('select goID = a.accID, t._Object_key, ec.abbreviation, refID = r.accID ' + \
	'from VOC_Annot t, ACC_Accession a, VOC_Evidence e, VOC_Term ec, ACC_Accession r ' + \
	'where t._AnnotType_key = 1000 ' + \
	'and t._Term_key = a._Object_key ' + \
	'and a._MGIType_key = 13 ' + \
	'and a.preferred = 1 ' + \
	'and t._Annot_key = e._Annot_key ' + \
	'and e._EvidenceTerm_key = ec._Term_key ' + \
	'and e._Refs_key in (61933,73197,73199) ' + \
	'and e._Refs_key = r._Object_key ' + \
	'and r._MGIType_key = 1 ' + \
	'and r._LogicalDB_key = 1 ' + \
	'and r.prefixPart= "J:"', 'auto')
for r in results:

    key = r['_Object_key']

    if key not in annot:
	annot[key] = []
    annot[key].append((r['goID'], r['abbreviation'], r['refID']))

    if r['goID'] not in annotByGOID:
        annotByGOID.append(r['goID'])

    if r['refID'] not in annotByRef:
        annotByRef.append(r['refID'])

#
# existing pubmed->j: relationships
#

results = db.sql('select jnumID = a1.accID, pubmedID = "PMID:" + a2.accID ' + \
	'from ACC_Accession a1, ACC_Accession a2 ' + \
	'where a1._MGIType_key = 1 ' + \
	'and a1._LogicalDB_key = 1 ' + \
	'and a1.prefixPart = "J:" ' + \
	'and a1.preferred = 1 ' + \
	'and a1._Object_key = a2._Object_key ' + \
	'and a2._MGIType_key = 1 ' + \
	'and a2._LogicalDB_key = 29 ' + \
	'and a2.preferred = 1 ', 'auto')
for r in results:
    key = r['pubmedID']
    value = r['jnumID']
    pubmed[key] = value

#
# GOA annotations
#

inFile = open(inFileName, 'r')

for line in inFile.readlines():

    tokens = string.split(line[:-1], '\t')
    databaseID = tokens[0]
    goaID = tokens[1]		# translate to MGI value
    goaSymbol = tokens[2]	# translate to MGI value
    notValue = tokens[3]
    goID = tokens[4]
    refID = tokens[5]		# translate to MGI value
    checkrefID = refID
    evidence = tokens[6]
    inferredFrom = tokens[7]
    dag = tokens[8]
    goaName = tokens[9]		# translate to MGI value
    synonyms = tokens[10]
    markerType = tokens[11]	# translate to MGI value
    taxID = tokens[12]
    modDate = tokens[13]
    assignedBy = tokens[14]

    # skip it if it's assigned by MGI; that means it came from us to begin with

    if assignedBy == 'MGI':
	mgiErrorFile.write(line)
	continue

    #
    # translate GOA "Refs" to MGI J: so we can check for duplicates
    #

    if refID == 'GOA:interpro':
	checkrefID = 'J:72247'

    if refID == 'GOA:spkw':
	checkrefID = 'J:60000'

    if refID == 'GOA:spec':
	checkrefID = 'J:72245'

    # error if GOA id is not found in MGI

    if not assoc.has_key(goaID):
	unresolvedBErrorFile.write(line)
	continue
    else:
        # error if GOA id maps to more than one MGI Marker

        if len(assoc[goaID]) > 1:
	    unresolvedAErrorFile.write(line)
	    continue

        mgiID = assoc[goaID][0]

    m = marker[mgiID]
    markerKey = m['_Marker_key']

    # duplicate error if the annotation already exists in MGI

    if annot.has_key(markerKey) and goID in annotByGOID and checkrefID in annotByRef:
        goaAnnot = (goID, evidence, checkrefID)
        if goaAnnot in annot[markerKey]:
	    dupErrorFile.write(line)
	    continue

    mgiFile.write(mgiLine % (databaseID, m['mgiID'], m['symbol'], notValue, goID, refID, evidence, inferredFrom,\
	dag, m['name'], synonyms, m['markerType'], taxID, modDate, assignedBy))

    # resolve pubmed ID to MGI J:

    if not pubmed.has_key(refID):
        if string.find(refID, 'PMID:') >= 0:
	    if evidence == 'ISS':
	        pubmedISSErrorFile.write(line)
            else:
	        pubmedErrorFile.write(line)
	continue
    else:
	jnumID = pubmed[refID]

    notValue = re.sub('_', ' ', notValue)
    annotFile.write(annotLine % (goID, m['mgiID'], jnumID, evidence, notValue, createdBy, loaddate))

inFile.close()

reportlib.finish_nonps(unresolvedAErrorFile)
reportlib.finish_nonps(unresolvedBErrorFile)
reportlib.finish_nonps(mgiErrorFile)
reportlib.finish_nonps(pubmedErrorFile)
reportlib.finish_nonps(pubmedISSErrorFile)
reportlib.finish_nonps(dupErrorFile)
reportlib.finish_nonps(mgiFile)
reportlib.finish_nonps(annotFile)

