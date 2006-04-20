#!/usr/local/bin/python

'''
#
# Purpose:
#
#	To load Annotations for the specified Vocabulary and MGI Object
#	into the VOC Annotation structures:  VOC_Annot, VOC_Evidence
#
#	See TR 2867.
#  
# Assumes:
#
#	That the Annotation Type (VOC_AnnotType) record already exists 
#	for the specified Vocabulary and MGI Object.
#
#	That no one else is adding Annotations to the database.
#
# Side Effects:
#
#	None
#
# Input:
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
#		field 10: Logical DB Name of Object (field 2), being Annotated (default is MGI)
#
# Parameters:
#	-S = database server
#	-D = database
#	-U = user
#	-P = password file
#	-M = mode (new, append, preview)
#	-I = input file
#	-A = annotation type name (ex. "PhenoSlim/Genotype", "GO/Marker")
#             the annotation type implicity defines the term vocabulary, 
#             the evidence vocabulary and the qualifier vocabulary
#	-R = reference (J:####) (used when mode = "new")
#	-O = load annotations to obsolete terms (default is to NOT load them)
#
#	processing modes:
#		new - delete the Annotations for the given Reference and Annotation Type.
#		      if reference = 0, then delete all Annotations
#		      of the specified Annotation Type
#
#		append - add Annotations from the input file (does not check for duplicates)
#
#		preview - perform all record verifications but do not load the data or
#		          make any changes to the database.  used for testing or to preview
#			  the load.
#
# Output:
#
#	1.  BCP files for VOC_Annot and VOC_Evidence
#	2.  Diagnostics file of all input parameters and SQL commands
#	3.  Error file
#
# Processing:
#
#	1. Verify the Annotation Type Name entry exists in VOC_AnnotType.
#	   If it does not exist, abort the load.
#	   If it does exist, retrieve the _AnnotType_key.
#
#	2. Verify Mode.
#		if mode = new and reference > 0: delete all Annotations for Annotation Type/Reference.
#		if mode = new and reference = 0: delete all Annotations for Annotation Type.
#		if mode = append: do nothing special
#		if mode = preview:  set "DEBUG" to True
#
#	3. Load Evidence Codes and Terms into dictionaries for quicker lookup.
#	    If not -O, then only load non-obsolete Terms.
#
#	For each line in the input file:
#
#	1.  Verify the Term Accession ID is valid for the given Annotation Type.
#	    If the verification fails, report the error and skip the record.
#
#	2.  Verify the MGI ID is valid for the given Annotation Type.
#	    If the verification fails, report the error and skip the record.
#	    If the verification succeeeds, store the MGI ID/Key pair in a dictionary
#	    for future reference.
#
#	3.  Verify the J: is valid.
#	    If the verification fails, report the error and skip the record.
#	    If the verification succeeeds, store the Jnum/Key pair in a dictionary
#	    for future reference.
#
#	4.  Verify the Evidence Code is valid for the given Annotation Type.
#	    If the verification fails, report the error and skip the record.
#
#	5.  Verify the Editor is provided (i.e. is not null).
#	    If the verification fails, report the error and skip the record.
#
#	6.  If the Date is not given, use the current date.
#
#	7.  Determine if an Annotation record for the given Term, Object and "Not" value
#	    already exists in VOC_Annot.  If it does exist, use its _Annot_key for the
#	    Evidence record.  If it does not exist, create a new Annotation record.
#
#	8.  Determine if an Evidence record for the given Annotation key, Evidence Code
#	    and Reference already exists in VOC_Evidence. If it does exist, write
#	    a message to the error file.  If it does not exist, create a new
#	    Evidence record.
#
# History:
#
# lec	10/04/2005
#	- TR 5188; replace isNot with _Qualifier_key in VOC_Annot.
#
# lec	01/28/2004
#	- TR 3404/JSAM
#
# lec	06/14/2002
#	- TR 3744 - added -O option to allow load of annotations obsolete terms
#		(the default is to allow load of annotations to obsolete terms)
#
# lec	05/28/2002
#	- TR 3724; fix deletion (add deletion of orphan VOC_Annot records)
#
# lec	01/22/2002
#	- created
#
'''

import sys
import os
import string
import getopt
import re
import db
import accessionlib
import mgi_utils
import loadlib
import vocabloadlib

#globals

DEBUG = 0		# set DEBUG to false unless preview mode is selected

inputFile = ''		# file descriptor
diagFile = ''		# file descriptor
errorFile = ''		# file descriptor
annotFile = ''		# file descriptor
evidenceFile = ''	# file descriptor
noteFile = ''		# file descriptor
noteChunkFile = ''	# file descriptor

diagFileName = ''	# file name
errorFileName = ''	# file name
annotFileName = ''	# file name
evidenceFileName = ''	# file name
noteFileName = ''	# file name
noteChunkFileName=  ''	# file name
passwordFileName = ''	# file name

mode = ''		# processing mode
delReference = 0	# deletion reference (J:###)
loadObsolete = 0	# load annotations to obsolete terms?
annotTypeName = ''	# VOC_AnnotType.name
annotTypeKey = 0	# VOC_AnnotType._AnnotType_key
annotKey = 0		# VOC_Annot._Annot_key
evidencePrimaryKey = 0	# VOC_Evidence._AnnotEvidence_key
noteKey = 0		# MGI_Note._Note_key
logicalDBKey = 1	# ACC_Accession._LogicalDB_key (default is "MGI", 1)
mgiNoteObjectKey = 25	# MGI_Note._MGIType_key
mgiNoteTypeKey = 1008	# MGI_Note._NoteType_key
mgiNoteSeqNum = 1	# MGI_NoteChunk.sequenceNum

termDict = {}		# dictionary of terms for quick lookup
objectDict = {}		# dictionary of objects for quick lookup
referenceDict = {}	# dictionary of references for quick lookup
ecodesDict = {}		# dictionary of evidence codes for quick lookup
annotDict = {}		# dictionary of annotation records for quick lookup
evidenceDict = {}	# dictionary of evidence records for quick lookup

loaddate = loadlib.loaddate

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
		'-M mode\n' + \
		'-I input file\n' + \
		'-A annotation type name\n' + \
		'-R reference\n' + \
		'-O load annotations to obsolete terms\n'
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
		diagFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
		errorFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
		diagFile.close()
		errorFile.close()
	except:
		pass

	db.useOneConnection()
	sys.exit(status)
 
def init():
	'''
	# requires: 
	#
	# effects: 
	# 1. Processes command line options
	# 2. Initializes local DBMS parameters
	# 3. Initializes global file descriptors/file names
	# 4. Initializes global keys
	#
	# returns:
	#
	'''
 
	global inputFile, diagFile, errorFile, errorFileName, diagFileName
	global annotFile, annotFileName, evidenceFile, evidenceFileName, passwordFileName
	global noteFile, noteFileName, noteChunkFile, noteChunkFileName
	global delReference, loadObsolete, mode
	global annotTypeKey, annotKey, annotTypeName, evidencePrimaryKey, noteKey
 
	try:
		optlist, args = getopt.getopt(sys.argv[1:], 'S:D:U:P:M:I:A:R:O')
	except:
		showUsage()
 
	#
	# Set server, database, user, passwords depending on options
	# specified by user.
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
			passwordFileName = opt[1]
                elif opt[0] == '-M':
                        mode = opt[1]
                elif opt[0] == '-I':
                        inputFileName = opt[1]
                elif opt[0] == '-A':
                        annotTypeName = opt[1]
                elif opt[0] == '-R':
                        delReference = opt[1]
                elif opt[0] == '-O':
                        loadObsolete = 1
                else:
                        showUsage()
 
	# Initialize db.py DBMS parameters
        password = string.strip(open(passwordFileName, 'r').readline())
	db.set_sqlLogin(user, password, server, database)
	db.useOneConnection(1)
 
	fdate = mgi_utils.date('%m%d%Y')	# current date
	head, tail = os.path.split(inputFileName) 
	diagFileName = tail + '.' + fdate + '.diagnostics'
	errorFileName = tail + '.' + fdate + '.error'
	annotFileName = tail + '.VOC_Annot.bcp'
	evidenceFileName = tail + '.VOC_Evidence.bcp'
	noteFileName = tail + '.MGI_Note.bcp'
	noteChunkFileName = tail + '.MGI_NoteChunk.bcp'

	try:
		inputFile = open(inputFileName, 'r')
	except:
		exit(1, 'Could not open file %s\n' % inputFileName)
		
	try:
		diagFile = open(diagFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % diagFileName)
		
	try:
		errorFile = open(errorFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % errorFileName)
		
	try:
		annotFile = open(annotFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % annotFileName)
		
	try:
		evidenceFile = open(evidenceFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % evidenceFileName)
		
	try:
		noteFile = open(noteFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % noteFileName)
		
	try:
		noteChunkFile = open(noteChunkFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % noteChunkFileName)
		
	# Log all SQL
	db.set_sqlLogFunction(db.sqlLogAll)

	# Set Log File Descriptor
	db.set_sqlLogFD(diagFile)

	diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
	diagFile.write('Server: %s\n' % (server))
	diagFile.write('Database: %s\n' % (database))
	diagFile.write('User: %s\n' % (user))
	diagFile.write('Annotation Type Name: %s\n' % (annotTypeName))
	diagFile.write('Annotation File: %s\n' % (inputFileName))
	diagFile.write('Deletion Reference: %s\n\n' % (delReference))

	errorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

def verifyAnnotType():
	'''
	# requires:
	#
	# effects:
	#	verifies that the Annotation Type Name exists in the VOC_AnnotType table
	#	if it does not exist, an error is written to the error file and the
	#	program is aborted.
	#	if it does exist, the global annotTypeKey is set to the primary key
	#
	# returns:
	#	nothing
	#
	'''

	global annotTypeKey, annotTypeName

	annotTypeName = re.sub('"', '', annotTypeName)

	results = db.sql('select _AnnotType_key from VOC_AnnotType where name = "%s"' % (annotTypeName), 'auto')

	if len(results) == 0:
		exit(1, 'Invalid Annotation Type Name: %s\n' % (annotTypeName))
	else:
		annotTypeKey = results[0]['_AnnotType_key']

def verifyMode():
	'''
	# requires:
	#
	# effects:
	#	Verifies the processing mode is valid.  If it is not valid,
	#	the program is aborted.
	#	Sets globals based on processing mode.
	#	Deletes data based on processing mode.
	#
	# returns:
	#	nothing
	#
	'''

	global DEBUG

	if mode == 'new':

		# verify deletion reference

		if delReference != "J:0":
			referenceKey = loadlib.verifyReference(delReference, 0, errorFile)

			if referenceKey is None:
				exit(1, 'Invalid Reference: %s\n' % (delReference))
		
			db.sql('delete VOC_Evidence from VOC_Annot a, VOC_Evidence e ' + \
				'where e._Refs_key = %s ' % (referenceKey) + \
				'and e._Annot_key = a._Annot_key ' + \
				'and a._AnnotType_key = %s\n' % (annotTypeKey), None, execute = not DEBUG)
			db.sql('delete VOC_Annot from VOC_Annot a ' + \
				'where a._AnnotType_key = %s ' % (annotTypeKey) + \
				'and not exists (select 1 from VOC_Evidence e ' + \
				'where a._Annot_key = e._Annot_key)', None, execute = not DEBUG)
		else:
			db.sql('delete VOC_Annot from VOC_Annot ' + \
				'where _AnnotType_key = %s\n' % (annotTypeKey), None, execute = not DEBUG)

	elif mode == 'append':
		pass

	elif mode == 'preview':
		DEBUG = 1
	else:
		exit(1, 'Invalid Processing Mode:  %s\n' % (mode))

def verifyTerm(termID, lineNum):
	'''
	# requires:
	#	termID - the Accession ID of the term
	#	lineNum - the line number of the record from the input file
	#
	# effects:
	#	verifies that:
	#		the Term exists 
	#		is of the appropriate type for the Annotation Type 
	#			of the load by checking the termDict dictionary for the term ID.
	#	writes to the error file if the Term is invalid
	#
	# returns:
	#	0 if the Term is invalid
	#	Term Key if the Term is valid
	#
	'''

	termKey = 0

	if termDict.has_key(termID):
		termKey = termDict[termID]
	else:
		if not loadObsolete:
			errorFile.write('Invalid or Obsolete Term (%d) %s\n' % (lineNum, termID))
		else:
			errorFile.write('Invalid Term (%d) %s\n' % (lineNum, termID))
		termKey = 0

	return(termKey)

def verifyObject(objectID, logicalDBKey, lineNum):
	'''
	# requires:
	#	objectID - the ID of the Object
	#	logicalDBKey - the logical DB of the objectID
	#	lineNum - the line number of the record from the input file
	#
	# effects:
	#	verifies that the Object exists and is of the appropriate type
	#	for the Annotation Type of the load by checking the objectDict
	#	dictionary for the object ID or the database.
	#	writes to the error file if the Object is invalid
	#	adds the Object ID/Key to the global objectDict dictionary if the
	#	object is valid
	#
	# returns:
	#	0 if the Object is invalid
	#	Object Key if the Object is valid
	#
	'''

	global objectDict

	if objectDict.has_key(objectID):
		objectKey = objectDict[objectID]
	else:
		results = db.sql('select a._Object_key ' + \
			'from ACC_Accession a, VOC_AnnotType t ' + \
			'where a.accID = "%s" ' % (objectID) + \
			'and a._LogicalDB_key = %s ' % (logicalDBKey) + \
			'and a._MGIType_key = t._MGIType_key ' + \
			'and t._AnnotType_key = %s\n' % (annotTypeKey), 'auto')

		if len(results) == 1:
			objectKey = results[0]['_Object_key']
			objectDict[objectID] = objectKey
		else:
			errorFile.write('Invalid Object (%d) %s\n' % (lineNum, objectID))
			objectKey = 0

	return(objectKey)

def setPrimaryKeys():
	'''
	# requires:
	#
	# effects:
	#	Sets the global primary keys values needed for the load
	#
	# returns:
	#	nothing
	#
	'''

	global annotKey, evidencePrimaryKey, noteKey

        results = db.sql('select maxKey = max(_Annot_key) + 1 from VOC_Annot', 'auto')
        if results[0]['maxKey'] is None:
                annotKey = 1000
        else:
                annotKey = results[0]['maxKey']

        results = db.sql('select maxKey = max(_AnnotEvidence_key) + 1 from VOC_Evidence', 'auto')
        if results[0]['maxKey'] is None:
                evidencePrimaryKey = 1000
        else:
                evidencePrimaryKey = results[0]['maxKey']

        results = db.sql('select maxKey = max(_Note_key) + 1 from MGI_Note', 'auto')
        if results[0]['maxKey'] is None:
                noteKey = 1000
        else:
                noteKey = results[0]['maxKey']

def loadDictionaries():
	'''
	# requires:
	#
	# effects:
	#	loads global dictionaries: ecodesDict, termDict
	#	for quicker lookup
	#
	# returns:
	#	nothing
	'''

	global ecodesDict, termDict

	results = db.sql('select e._Term_key, e.abbreviation ' + \
			'from VOC_Term e, VOC_AnnotType t ' + \
			'where e._Vocab_key = t._EvidenceVocab_key ' + \
			'and t._AnnotType_key = %s\n' % (annotTypeKey), 'auto')

	for r in results:
		ecodesDict[r['abbreviation']] = r['_Term_key']

	cmd = 'select t._Object_key, t.accID ' + \
		'from VOC_Term_Acc_View t, VOC_Term tm, VOC_AnnotType a ' + \
		'where t._Object_key = tm._Term_key ' + \
		'and tm._Vocab_key = a._Vocab_key ' + \
		'and a._AnnotType_key = %s\n' % (annotTypeKey)

	# if loadObsolete is false, then only load non-obsoleted terms...

	if not loadObsolete:
		# only load non-obsoleted terms
		cmd = cmd + 'and tm.isObsolete = 0'

	results = db.sql(cmd, 'auto')

	for r in results:
		termDict[r['accID']] = r['_Object_key']

def createAnnotationRecord(objectKey, termKey, qualifierKey, entryDate):
	'''
	# requires:
	#	objectKey - primary key of the Object
	#	termKey - primary key of the Term
	#	qualifierKey - primary key of the Qualifier
	#	entryDate - creation and modification date of Annotation
	#
	# effects:
	#	determines the primary key of the Annotation record
	#	by checking the global annotDict dictionary or the database
	#
	#	if this is a new Annotation record, writes a new record
	#	to the annotFile (bcp), adds the new annotation key
	#	to the global annotDict dictionary and increments the global
	#	annotKey counter.
	#
	# returns:
	#	primary key for Annotation record
	#
	'''

	global annotKey, annotDict

	# if an annotation already exists for the same Object/Term/Not, 
	# use the same annotation key

	aKey = '%s:%s:%s:%s' % (annotTypeKey, objectKey, termKey, qualifierKey)

	# annotation may exist in our dictionary already...

	if annotDict.has_key(aKey):
		useAnnotKey = annotDict[aKey]
	else:
		# not in our dictionary, let's try to find it in the database

		results = db.sql('select _Annot_key from VOC_Annot ' + \
			'where _AnnotType_key = %s ' % (annotTypeKey) + \
			'and _Object_key = %s ' % (objectKey) + \
			'and _Term_key = %s ' % (termKey) + \
			'and _Qualifier_key = %s ' % (qualifierKey), 'auto')
			
		# found it in the database

		if len(results) > 0:
			useAnnotKey = results[0]['_Annot_key']
			annotDict[aKey] = useAnnotKey

		# not found in the database; let's create it

		else:
			useAnnotKey = annotKey
			annotDict[aKey] = useAnnotKey
			annotKey = annotKey + 1

			# create the new VOC_Annot record

			annotFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
			% (useAnnotKey, annotTypeKey, objectKey, termKey, \
		   	qualifierKey, entryDate, entryDate))

	return(useAnnotKey)

def createEvidenceRecord(newAnnotKey, evidenceKey, referenceKey, inferredFrom, editorKey, notes, entryDate, lineNum):
	'''
	# requires:
	#	newAnnotKey - primary key of the Annotation object
	#	evidenceKey - primary key of the Evidence Code
	#	referenceKey - primary key of the Reference
	#	inferredFrom - inferred from value
	#	editorKey - primary key of the Editor
	#	notes - notes
	#	entryDate - creation and modification date of Annotation
	#	lineNum - the line number of the record from the input file
	#
	# effects:
	#	determines if this Evidence record is a duplicate
	#	by checking the global evidenceDict dictionary or the database
	#
	#	if this is a duplicate, then writes to the error file and returns
	#	if new Evidence record, writes a new record to the evidenceFile (bcp)
	#	and adds the new evidence key to the global evidenceDict dictionary
	#
	# returns:
	#	nothing
	#
	'''

	global evidencePrimaryKey, evidenceDict, noteKey

	# make sure this is not a duplicate evidence statement

	eKey = '%s:%s:%s' % (newAnnotKey, evidenceKey, referenceKey)

	# evidence record may exist in our dictionary already
	# if so, it's a duplicate; let's report it

	if evidenceDict.has_key(eKey):
		errorFile.write('Duplicate Evidence Statement (in input file): %d\n' % (lineNum))
		return

	# not a duplicate

	evidenceDict[eKey] = eKey

	# not in our dictionary, let's try to find it in the database

	results = db.sql('select _Annot_key from VOC_Evidence ' + \
		'where _Annot_key = %s ' % (newAnnotKey) + \
		'and _EvidenceTerm_key = %s ' % (evidenceKey) + \
		'and _Refs_key = %s ' % (referenceKey), 'auto')
			
	# found it in the database; it's a duplicate

	if len(results) > 0:
		errorFile.write('Duplicate Evidence Statement (in database already): %d\n' % (lineNum))
		return

	# not found in the database; let's create it

	evidenceFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
		% (evidencePrimaryKey, newAnnotKey, evidenceKey, referenceKey, inferredFrom, \
		   editorKey, editorKey, entryDate, entryDate))

	mgiNoteSeqNum = 1
	if len(notes) > 0:
	    noteFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
		% (noteKey, evidencePrimaryKey, mgiNoteObjectKey, mgiNoteTypeKey, editorKey, editorKey, entryDate, entryDate))
	    while len(notes) > 255:
	        noteChunkFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
		    % (noteKey, mgiNoteSeqNum, notes[:255], editorKey, editorKey, entryDate, entryDate))
		notes = notes[255:]
		mgiNoteSeqNum = mgiNoteSeqNum + 1

	    if len(notes) > 0:
	        noteChunkFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
		    % (noteKey, mgiNoteSeqNum, notes, editorKey, editorKey, entryDate, entryDate))
	    noteKey = noteKey + 1

	evidencePrimaryKey = evidencePrimaryKey + 1

def processFile():
	'''
	# requires:
	#
	# effects:
	#	Reads input file
	#	Verifies and Processes each line in the input file
	#
	# returns:
	#	nothing
	#
	'''

	global logicalDBKey

	lineNum = 0

	# For each line in the input file

	for line in inputFile.readlines():

		error = 0
		lineNum = lineNum + 1

		# Split the line into tokens
		tokens = string.splitfields(line[:-1], '\t')

		try:
			termID = tokens[0]
			objectID = tokens[1]
			jnum = tokens[2]
			evidence = tokens[3]
			inferredFrom = string.strip(tokens[4])
			qualifier = string.strip(tokens[5])
			editor = string.strip(tokens[6])
			entryDate = string.strip(tokens[7])
			notes = string.strip(tokens[8])

			if len(tokens) == 10:
				logicalDBKey = accessionlib.get_LogicalDB_key(tokens[9])

		except:
			exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

		termKey = verifyTerm(termID, lineNum)
		objectKey = verifyObject(objectID, logicalDBKey, lineNum)
		referenceKey = loadlib.verifyReference(jnum, lineNum, errorFile)
		evidenceKey = vocabloadlib.verifyEvidence(evidence, annotTypeKey, lineNum, errorFile)
		qualifierKey = vocabloadlib.verifyQualifier(qualifier, annotTypeKey, lineNum, errorFile)
		editorKey = loadlib.verifyUser(editor, lineNum, errorFile)

		if termKey == 0 or objectKey == 0 or \
			referenceKey == 0 or \
			evidenceKey == 0 or \
			qualifierKey == 0 or \
			editorKey == 0:

			# set error flag to true
			error = 1

		# If the entry date is not given, use the current date

		if len(entryDate) == 0:
			entryDate = loaddate

		# if errors, continue to next record
		if error:
			continue

		# if no errors, process the annotation

		newAnnotKey = createAnnotationRecord(objectKey, termKey, qualifierKey, entryDate)

		createEvidenceRecord(newAnnotKey, \
			evidenceKey, \
			referenceKey, \
			inferredFrom, \
			editorKey, \
			notes, \
			entryDate, \
			lineNum)

#	end of "for line in inputFile.readlines():"

def bcpFiles():
	'''
	# requires:
	#
	# effects:
	#	BCPs the data into the database
	#
	# returns:
	#	nothing
	#
	'''

	if DEBUG:
		return

	annotFile.close()
	bcpAnnot = 'cat %s | bcp %s..%s in %s -c -t\"\t" -e %s -S%s -U%s >> %s' \
		% (passwordFileName, db.get_sqlDatabase(), \
	   	'VOC_Annot', annotFileName, errorFileName, db.get_sqlServer(), db.get_sqlUser(), diagFileName)
	diagFile.write('%s\n' % bcpAnnot)

	evidenceFile.close()
	bcpEvidence = 'cat %s | bcp %s..%s in %s -c -t\"\t" -e %s -S%s -U%s >> %s' \
		% (passwordFileName, db.get_sqlDatabase(), \
	   	'VOC_Evidence', evidenceFileName, errorFileName, db.get_sqlServer(), db.get_sqlUser(), diagFileName)
	diagFile.write('%s\n' % bcpEvidence)

	noteFile.close()
	bcpNote = 'cat %s | bcp %s..%s in %s -c -t\"\t" -e %s -S%s -U%s >> %s' \
		% (passwordFileName, db.get_sqlDatabase(), \
	   	'MGI_Note', noteFileName, errorFileName, db.get_sqlServer(), db.get_sqlUser(), diagFileName)
	diagFile.write('%s\n' % bcpNote)

	noteChunkFile.close()
	bcpNoteChunk = 'cat %s | bcp %s..%s in %s -c -t\"\t" -e %s -S%s -U%s >> %s' \
		% (passwordFileName, db.get_sqlDatabase(), \
	   	'MGI_NoteChunk', noteChunkFileName, errorFileName, db.get_sqlServer(), db.get_sqlUser(), diagFileName)
	diagFile.write('%s\n' % bcpNoteChunk)

	os.system(bcpAnnot)
	os.system(bcpEvidence)
	os.system(bcpNote)
	os.system(bcpNoteChunk)

#	db.sql('dump transaction %s with truncate_only' % (db.get_sqlDatabase()), None)

#
# Main
#

init()
verifyAnnotType()
verifyMode()
setPrimaryKeys()
loadDictionaries()
processFile()
bcpFiles()
exit(0)

