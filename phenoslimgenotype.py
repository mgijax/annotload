#!/usr/local/bin/python

'''
#
# Purpose:
#
#	To translate a PhenoSlim Annotation Spreadsheet into an input
#	file for the Annotation Load.
#
#	See TRs 2239 and 2867.
#  
# 	Used only for migration of Spreadsheets into database.
#	Used for MGI 3.8 migration only.
#
# Assumes:
#
#	All Genotypes (Strain/Gene/Allele 1/Allele 2) are single mutants.
#
# Side Effects:
#
#	Creates Strain records for new Strains.
#	Creates Genotype records for new Strain/Gene/Allele 1/Allele 2 combinations.
#
# Input:
#
#	A tab-delimited file in the format:
#		field 1: MGI ID of Gene
#		field 2: Gene Symbol
#		field 3: MGI ID of Allele 1
#		field 4: Allele 1 Symbol
#		field 5: MGI ID of Allele 2
#		field 6: Allele 2 Symbol
#		field 7: Strain
#		field 8: PhenoSlim Term
#		field 9: J:
#		field 10: Evidence Code
#		field 11: Notes
#
# Parameters:
#	-S = database server
#	-D = database
#	-U = user
#	-P = password file
#	-I = input file
#	-E = editor
#
# Output:
#
#	1.  Tab-delimited file in the format:
#		field 1: Accession ID of PhenoSlim Term
#		field 2: MGI ID of Genotype
#		field 3: J:
#		field 4: Evidence Code Abbreviation
#		field 5: Inferred From
#		field 6: NOT
#		field 7: Editor
#		field 8: Date
#		field 9: Notes
#
# Processing:
#
#	For each line in the input file:
#
#	1.  Retrieve the primary key of the Marker by MGI ID.
#
#	2.  Retrieve the primary key of each Allele by MGI ID.
#
#	3.  Retrieve the primary key of the Strain.
#
#	4.  Create the Genotype record (if it does not already exist) for the 
#	    Strain/Marker/Allele 1/Allele 2.
#
#	5.  Write the record to the output file.
#
# History:
#
# lec	01/22/2002
#	- created
#
'''

import sys
import os
import string
import getopt
import db
import mgi_utils

#globals

inputFile = ''
outputFile = ''
errorFile = ''

editor = ''
genotypeKey = 0
allelePairKey = 0
newStrainKey = 0

termDict = {}		# dictionary of ps term to speed up lookup

def showUsage():
	'''
	# requires:
	#
	# effects:
	# Displays the correct usage of this program and exits
	# with status of 1.
	#
	# returns:
	'''
 
	usage = 'usage: %s -S server\n' % sys.argv[0] + \
		'-D database\n' + \
		'-U user\n' + \
		'-P password file\n' + \
		'-I input file\n' + \
		'-E editor\n'
	exit(1, usage)
 
def exit(status, message = None):
	'''
	# requires: status, the numeric exit status (integer)
	#           message (string)
	#
	# effects:
	# Print message to stderr and exits
	#
	# returns:
	#
	'''
 
	if message is not None:
		sys.stderr.write('\n' + str(message) + '\n')
 
	try:
		outputFile.close()
		errorFile.close()
	except:
		pass

	sys.exit(status)
 
def init():
	'''
	# requires: 
	#
	# effects: 
	# 1. Processes command line options
	# 2. Initializes local DBMS parameters
	# 3. Initializes global file descriptors
	#
	# returns:
	#
	'''
 
	global inputFile, outputFile, errorFile
	global genotypeKey, allelePairKey, newStrainKey
	global editor
	global termDict
 
	try:
		optlist, args = getopt.getopt(sys.argv[1:], 'S:D:U:P:I:E:')
	except:
		showUsage()
 
	#
	# Set server, database depending on options specified by user.
	#
 
	server = None
	database = None
	user = None
	password = None
 
	for opt in optlist:
                if opt[0] == '-S':
                        server = opt[1]
                elif opt[0] == '-D':
                        database = opt[1]
                elif opt[0] == '-U':
                        user = opt[1]
                elif opt[0] == '-P':
                        password = string.strip(open(opt[1], 'r').readline())
                elif opt[0] == '-I':
                        inputFileName = opt[1]
                elif opt[0] == '-E':
                        editor = opt[1]
                else:
                        showUsage()
 
	# User must specify Server, Database, User and Password
	if server is None or database is None or user is None or password is None:
		showUsage()
 
	# Initialize db.py DBMS parameters
	db.set_sqlLogin(user, password, server, database)
 
	head, tail = os.path.split(inputFileName)
	outputFileName = tail + '.annotload'
	errorFileName = outputFileName + '.error'

	try:
		inputFile = open(inputFileName, 'r')
	except:
		exit(1, 'Could not open file %s\n' % inputFileName)
		
	try:
		outputFile = open(outputFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % outputFileName)
		
	try:
		errorFile = open(errorFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % errorFileName)
		
	# Get next available primary keys
	#

	results = db.sql('select maxKey = max(_Genotype_key) from GXD_Genotype', 'auto')
	if results[0]['maxKey'] is None:
		genotypeKey = 1000
	else:
		genotypeKey = results[0]['maxKey']

	results = db.sql('select maxKey = max(_AllelePair_key) from GXD_AllelePair', 'auto')
	if results[0]['maxKey'] is None:
		allelePairKey = 1000
	else:
		allelePairKey = results[0]['maxKey']

	results = db.sql('select maxKey = max(_Strain_key) from PRB_Strain', 'auto')
	if results[0]['maxKey'] is None:
		newStrainKey = 1000
	else:
		newStrainKey = results[0]['maxKey']

	results = db.sql('select accID, term from VOC_Term_View where _Vocab_key = 1', 'auto')
	for r in results:
		termDict[r['term']] = r['accID']

def processFile():
	'''
	# requires:
	#
	# effects:
	#	Creates output file
	#
	# returns:
	#	nothing
	#
	'''

	global genotypeKey, allelePairKey, newStrainKey

	# use dictionaries to store accID/key pairs for quicker lookup

	markerDict = {}
	alleleDict = {}
	strainDict = {}
	genotypeDict = {}

	# For each line in the input file
	lineNum = 0

	for line in inputFile.readlines():

		# Split the line into tokens

		lineNum = lineNum + 1
		error = 0
		tokens = string.splitfields(line[:-1], '\t')

		try:
			markerID = string.strip(tokens[0])
			marker = tokens[1]
			allele1ID = string.strip(tokens[2])
			allele1 = tokens[3]
			allele2ID = string.strip(tokens[4])
			allele2 = tokens[5]
			strain = string.strip(tokens[6])
			psterm = string.strip(tokens[7])
			jnum = string.strip(tokens[8])
			ecode = string.strip(tokens[9])
	
		except:
			errorFile.write('Problem with line %d.\n' % (lineNum))
			errorFile.write(str(tokens) + '\n\n')
			continue

		# for some reason, this field is not always present in the converted spreadsheet
		# some weird Excel thing i suppose...

		if len(tokens) == 11:
			if len(tokens[10]) > 255:
				notes = tokens[10][:254]
			else:
				notes = tokens[10]
		else:
			notes = ''

		if string.find(markerID, 'MGI:') < 0:
			markerID = 'MGI:' + markerID

		if string.find(allele1ID, 'MGI:') < 0:
			allele1ID = 'MGI:' + allele1ID

		if len(allele2ID) > 0 and string.find(allele2ID, 'MGI:') < 0:
			allele2ID = 'MGI:' + allele2ID

		if allele2ID = '?':
			isUnknown = 1
			allele2ID = ''
		else:
			isUnknown = 0

		# get marker key

		if markerDict.has_key(markerID):
			markerKey = markerDict[markerID]
		else:
			results = db.sql('select _Object_key from MRK_Acc_View where accid = "%s"' % (markerID), 'auto')
			if len(results) > 0:
				markerKey = results[0]['_Object_key']
				markerDict[markerID] = markerKey
			else:
				errorFile.write('Invalid Marker (%d): %s\n' % (lineNum, markerID))
				error = 1

		# get allele 1 key

		if alleleDict.has_key(allele1ID):
			allele1Key = alleleDict[allele1ID]
		else:
			results = db.sql('select _Object_key from ALL_Acc_View where accid = "%s"' % (allele1ID), 'auto')
			if len(results) > 0:
				allele1Key = results[0]['_Object_key']
				alleleDict[allele1ID] = allele1Key
			else:
				errorFile.write('Invalid Allele (%d): %s\n' % (lineNum, allele1ID))
				error = 1

		# get allele 2 key, if allele 2 exists

		if len(allele2ID) > 0:
			if alleleDict.has_key(allele2ID):
				allele2Key = alleleDict[allele2ID]
			else:
				results = db.sql('select _Object_key from ALL_Acc_View where accid = "%s"' % (allele2ID), 'auto')
				if len(results) > 0:
					allele2Key = results[0]['_Object_key']
					alleleDict[allele2ID] = allele2Key
				else:
					errorFile.write('Invalid Allele (%d): %s\n' % (lineNum, allele2ID))
					error = 1
		else:
			allele2Key = 'NULL'

		# get strain key; create a new strain if it doesn't exist

		if strainDict.has_key(strain):
			strainKey = strainDict[strain]
		else:
			results = db.sql('select _Strain_key from PRB_Strain where strain = "%s"' % (strain), 'auto')
			if len(results) > 0:
				strainKey = results[0]['_Strain_key']
				strainDict[strain] = strainKey
			else:
				# create the strain
				newStrainKey = newStrainKey + 1
				cmd = 'insert PRB_Strain values(%s, "%s", 0, 0, 0, getdate(), getdate())\n' % (newStrainKey, strain) + \
					'insert MLP_Strain values(%s, -1, NULL, NULL, getdate(), getdate())\n' % (newStrainKey) + \
					'insert MLP_Extra values(%s,NULL,NULL,NULL,NULL,getdate(),getdate())\n' % (newStrainKey)
				db.sql(cmd, None)
				strainKey = newStrainKey
				strainDict[strain] = strainKey

		# get term ID

		if termDict.has_key(psterm):
			termID = termDict[psterm]
		else:
			errorFile.write('Invalid Term (%d): %s\n' % (lineNum, psterm))
			error = 1

		# if no errors, we can look for the genotype record

		if not error:
			
			gKey = '%s:%s:%s' % (strain, allele1ID, allele2ID)

			if genotypeDict.has_key(gKey):
				genotypeID = genotypeDict[gKey]

			# add a genotype record (if one does not already exist)
			# we're only adding single mutants

			else:
				results = db.sql('select a.accID from GXD_AllelePair_View v, GXD_Genotype_Acc_View a ' + \
					'where v._Marker_key = %s ' % (markerKey) + \
					'and v._Allele_key_1 = %s ' % (allele1Key) + \
					'and v._Allele_key_2 = %s ' % (allele2Key) + \
					'and v._Strain_key = %s ' % (strainKey) + \
					'and v._Genotype_key = a._Object_key ' + \
					'and not exists (select 1 from GXD_AllelePair v2 ' + \
					'where v._Genotype_key = v2._Genotype_key ' + \
					'and v2.sequenceNum > 1)', 'auto')

				# if genotype record does not exist

				if len(results) == 0:

					# increment the primary key counters
					genotypeKey = genotypeKey + 1
					allelePairKey = allelePairKey + 1

					# create the record
					cmd = 'insert into GXD_Genotype ' + \
						'values(%s, %s, 0, "%s", "%s", NULL, getdate(), getdate())\n' \
						% (genotypeKey, strainKey, editor, editor) + \
						'insert into GXD_AllelePair ' + \
						'values (%s, %s, 1, %s, %s, %s, %d, getdate(), getdate())\n' \
						% (allelePairKey, genotypeKey, allele1Key, allele2Key, markerKey, isUnknown)
					db.sql(cmd, None)

					# grab the MGI Acc ID of the new record
					results = db.sql('select accID from GXD_Genotype_Acc_View ' + \
						'where _Object_key = %s' % (genotypeKey), 'auto')
					if len(results) > 0:
						genotypeID = results[0]['accID']
						genotypeDict[gKey] = genotypeID
					else:
						errorFile.write('Could Not Create Genotype Record for Strain: %s\n') % (strain)
						error = 1

				# single mutant exists
				elif len(results) == 1:
					# genotype exists
					genotypeID = results[0]['accID']
					genotypeDict[gKey] = genotypeID

			# write output file

			if not error:

				outputFile.write('%s\t%s\t%s\t%s\t\t\t%s\t\t%s\n' \
					% (termID, genotypeID, jnum, ecode, editor, mgi_utils.prvalue(notes)))


#
# Main
#

init()
processFile()
exit(0)

