#!/usr/local/bin/python

'''
#
# Purpose:
#
#	To translate the GO Annotation Spreadsheet into an input
#	file for the Annotation Load.
#
#	See TR 2867.
#  
# Assumes:
#
# Side Effects:
#
# Input:
#
#	A tab-delimited file in the format:
#		field 1: Gene Symbol
#		field 2: Gene Name
#		field 3: MGI ID of Gene
#		field 4: Vocabulary DAG (ex. GO Function, GO Cellular Component)
#		field 5: Not
#		field 6: GO ID
#		field 7: GO Definition
#		field 8: Reference (J:#####)
#		field 9: Evidence
#		field 10: Inferred From
#		field 11: Editor
#		field 12: Date
#		field 13: Notes
#		field 14: Confirmed
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
#	1.  Retrieve Evidence Code abbreviation for given Evidence
#
#	2.  Write the record to the output file.
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
 
	try:
		optlist, args = getopt.getopt(sys.argv[1:], 'S:D:I:O:')
	except:
		showUsage()
 
	#
	# Set server, database, depending on options specified by user.
	#
 
	server = None
	database = None
 
	for opt in optlist:
                if opt[0] == '-S':
                        server = opt[1]
                elif opt[0] == '-D':
                        database = opt[1]
                elif opt[0] == '-I':
                        inputFileName = opt[1]
                else:
                        showUsage()
 
	# User must specify Server, Database
	if server is None or database is None:
		showUsage()
 
	# Initialize db.py DBMS parameters
	db.set_sqlServer(server)
	db.set_sqlDatabase(database)
 
	outputFileName = inputFileName + ".annotload"
	errorFileName = outputFileName + ".error"

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

	# For each line in the input file

	for line in inputFile.readlines():

		# Split the line into tokens

		error = 0
		tokens = string.splitfields(line[:-1], '\t')

		# we're only assigning the tokens we're interested in

		markerID = string.strip(tokens[2])
		notTerm = string.strip(tokens[4])
		goID = string.strip(tokens[5])
		jnum = string.strip(tokens[7])
		evidence = string.strip(tokens[8])
		inferredFrom = regsub.gsub('|', ', ', string.strip(tokens[9]))
		editor = string.strip(tokens[10])
		entryDate = string.strip(tokens[11])
		notes = string.strip(tokens[12])

		# if term = unknown, then set Evidence Code = NAS and J = ???
		# and move original J: into Notes

		if goID in ["GO:0000004", "GO:0008372", "GO:0005554"]:
			evidence = "non-traceable author statement"
			notes = jnum
			jnum = "J:73796"

		# get evidence code

		if evidence == "author said so":
			ecode = "TAS"
		elif evidence == "non-traceable author statement":
			ecode = "NAS"
		elif evidence == "inferred from direct assay":
			ecode = "IDA"
		elif evidence == "inferred from electronic annotation":
			ecode = "IEA"
		elif string.find(evidence, "inferred from genetic interaction") > -1:
			ecode = "IGI"
		elif string.find(evidence, "inferred from physical interaction") > -1:
			ecode = "IPI"
		elif string.find(evidence, "inferred from physica? interaction") > -1:
			ecode = "IPI"
		elif evidence == "inferred from mutant phenotype":
			ecode = "IMP"
		elif evidence == "not available":
			ecode = "ND"
		elif string.find(evidence, "inferred from sequence similarity") > -1:
			ecode = "ISS"
		else:
			errorFile.write("Could Not Determine Evidence Code: " + evidence + "\n")
			error = 1

		if string.find(notes, 'Hand') > -1:
			i = string.find(notes, 'Hand')
			notes = notes[i + 5:-1]

		if notes in ['EC', 'GOFISH']:
			notes = None

		if notTerm == "Not":
			notTerm = "NOT"

		if entryDate == "Mar-00":
			entryDate = "03/01/00"

		# if no errors

		if not error:
			
			outputFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
				% (goID, markerID, jnum[2:], ecode, mgi_utils.prvalue(inferredFrom), \
					mgi_utils.prvalue(notTerm), editor, entryDate, mgi_utils.prvalue(notes)))

#
# Main
#

init()
processFile()
exit(0)

