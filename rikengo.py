#!/usr/local/bin/python

'''
#
# Purpose:
#
#	To translate the RIKEN/GO Annotations into an input file for the Annotation Load.
#	To generate a file of new genes (for rikenloadprepub)...
#
# Assumes:
#
# Side Effects:
#
# Input:
#
#	A tab-delimited file in the format:
#		field 1: Riken Clone ID
#		field 2: GO ID
#		field 3: |-delimited evidence codes 
#			the first one is the one we're interested in)
#		field 4: DB Reference
#			in GO-format 
#			(MGD|MGI:123345, SCOP|12345, SWISSPROT|Q9DD23, SPTR|P18572, InterPro|IPR000865)
#
# Parameters:
#	-S = database server
#	-D = database
#	-I = input file
#
# Output:
#
#	1.  Tab-delimited file in the format:
#		field 1: Accession ID of GO ID
#		field 2: MGI ID of Marker
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
#	1.  Find corresponding Marker MGI ID for given RIKEN Clone ID
#
#	2.  Extract the first evidence code (IEA or ND)
#
#	3.  Extract the DB reference and translate into MGI format:
#		(MGI:12345, SP:12345, etc., IP:IPR000865)
#
#	4.  Write the record to the output file.
#
# History:
#
# lec	07/01/2002
#	- created
#
'''

import sys
import os
import string
import getopt
import regsub
import db
import mgi_utils

#globals

inputFile = ''
outputFile = ''
errorFile = ''

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
		'-I input file\n'
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
		newgenesFile.close()
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
 
	global outputFile, errorFile, newgenesFile, inputFileName
 
	try:
		optlist, args = getopt.getopt(sys.argv[1:], 'S:D:U:P:I:O:')
	except:
		showUsage()
 
	#
	# Set server, database, depending on options specified by user.
	#
 
	server = None
	database = None
	user = None
	passwordFileName = None
 
	for opt in optlist:
                if opt[0] == '-S':
                        server = opt[1]
                elif opt[0] == '-D':
                        database = opt[1]
                elif opt[0] == '-U':
                        user = opt[1]
                elif opt[0] == '-P':
                        passwordFileName = opt[1]
                elif opt[0] == '-I':
                        inputFileName = opt[1]
                else:
                        showUsage()
 
	# User must specify Server, Database
	if server is None or database is None or user is None or passwordFileName is None:
		showUsage()
 
	# Initialize db.py DBMS parameters
	password = string.strip(open(passwordFileName, 'r').readline())
	db.set_sqlLogin(user, password, server, database)

	head, tail = os.path.split(inputFileName)
	outputFileName = tail + ".annotload"
	errorFileName = outputFileName + ".error"
	newgenesFileName = outputFileName + ".newgenes"

	try:
		outputFile = open(outputFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % outputFileName)
		
	try:
		errorFile = open(errorFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % errorFileName)
		
	try:
		newgenesFile = open(newgenesFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % newgenesFileName)
		
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

	# dictionary of clone id/marker id pairs
	cloneID = {}

	# get clone id/mgi ids from Fantom table or MGD table
	results = db.sql('select clone = riken_cloneID, mgiID = final_mgiID from MGI_Fantom2 ' + \
		'where riken_cloneID != "zilch" and final_mgiID like "MGI:%" ' + \
		'union ' + \
		'select clone = substring(symbol, 1, charindex("Rik", symbol) - 1), mgiID ' + \
		'from MRK_Mouse_View ' + \
		'where symbol like "%Rik"', 'auto')
	for r in results:
		cloneID[r['clone']] = r['mgiID']

	# get clone id/mgi ids from input file
	inputFile = open(inputFileName, 'r')
	for line in inputFile.readlines():
		tokens = string.splitfields(line[:-1], '\t')
		clone = string.strip(tokens[0])
		inferredFrom = regsub.gsub('MGD|', '', string.strip(tokens[3]))
		if not cloneID.has_key(clone):
			if string.find(inferredFrom, 'MGI:') > -1:
				cloneID[clone] = inferredFrom
	inputFile.close()

	newgenes = []
	entryDate = mgi_utils.date('%m/%d/%Y')
	editor = 'riken'
	jnum = '65060'

	inputFile = open(inputFileName, 'r')
	for line in inputFile.readlines():

		# Split the line into tokens

		error = 0
		tokens = string.splitfields(line[:-1], '\t')

		# we're only assigning the tokens we're interested in

		clone = string.strip(tokens[0])
		goID = string.strip(tokens[1])
		evidences = string.splitfields(string.strip(tokens[2]), '|')
		# we're only interested in the first evidence code
		ecode = evidences[0]

		if cloneID.has_key(clone):
			markerID = cloneID[clone]

		# add unique new "genes" to newgenes file for processing by rikenloadprepub
		else:
			if clone not in newgenes:
				newgenesFile.write('%s\t\t\t\n' % (clone))
				newgenes.append(clone)
			error = 1

		if not error:
			outputFile.write('%s\t%s\t%s\t%s\t\t\t%s\t%s\t\n' \
				% (goID, markerID, jnum, ecode, editor, entryDate))

#
# Main
#

init()
processFile()
exit(0)
