'''
#
# Purpose:
#
#    To load Annotations for the specified Vocabulary and MGI Object into the VOC Annotation structures:
#
#    VOC_Annot
#    VOC_Evidence
#    VOC_Evidence_Property (TR10044)
#    MGI_Note (_mgitype_key = 41 used by goload/bin/go_annot_extensions_display.sh go_isoforms_display.sh
#
# Assumes:
#
#    That the Annotation Type (VOC_AnnotType) record already exists for the specified Vocabulary and MGI Object.
#    That no one else is adding Annotations to the database.
#    Usage: annotload.py [mcv|mp|go|diseaseMarker|diseaseAllele|mpMarker|mpAllele|omimhpo]
#
# Input:
#
#    A tab-delimited file in the format:
#
#       field 1: Accession ID of Vocabulary Term being Annotated to
#       field 2: ID of MGI Object being Annotated (ex. MGI ID)
#       field 3: J: (J:#####)
#       field 4: Evidence Code Abbreviation (max length 5)
#       field 5: Inferred From 
#       field 6: Qualifier 
#       field 7: Editor (max length 30)
#       field 8: Date (MM/DD/YYYY)
#       field 9: Notes 
#
#       field 10: Logical DB Name of Object (field 2),
#               being Annotated (default is MGI); optional
#               field 10 is not required and should ONLY be used if
#               you need to specify a logoical DB key OTHER THAN "1" (MGI).
#               For example, see "entrezgeneload/human/annotations.py".
#
#       field 11: Properties
#           &=&   : separates property from its value 
#               prop1&=&val1
#
#           &==&  : separates one property/value pair from another
#               prop1&=&val1&==&prop2&=&val2
#
#           &===& : separates a stanza; stanza contain >= 1 property/value pair
#               prop1&=&val1&==&prop2&=&val2&===&prop1&=&val1&==&prop2&=&val2
#
#       field 11 GO examples:
#
#       one stanza, one property/value pair:
#           term&=&UniProtKB:P12023
#
#       one stanza, two property/value pair:
#           term&=&UniProtKB:P12023&==&external ref&=&PMID:2834384|EXP
#
#       two stanza, two property/value pair:
#           term&=&UniProtKB:P12023&==&external ref&=&PMID:2834384|EXP&===&term&=&UniProtKB:P12023&==&external ref&=&PMID:2834384|EXP
#
#	field 12+ may be used by curators and this load should ignore
#
# Parameters:
#
#    processing modes:
#       new - delete the Annots for the given Reference or User, and Annot Type.
#           if reference = 0 and user is none, then delete all Annots of the specified Annot Type
#
#       append - add Annots from the input file
#           (does not check for duplicates)
#
#       preview - perform all record verifications but do not load the data or
#           make any changes to the database. Used for testing or to preview the load.
#
#       delete - delete the Annots for the given Reference or User, and Annot Type.
#
# Output:
#
#    1.  BCP files for VOC_Annot, VOC_Evidence, VOC_Evidence_Property
#    2.  Diagnostics file of all input parameters and SQL commands
#    3.  Error file
#
# Processing:
#
#    1. Verify the Annotation Type Name entry exists in VOC_AnnotType.
#       If it does not exist, abort the load.
#       If it does exist, retrieve the _AnnotType_key.
#
#    2. Verify Mode.
#       if mode = new and reference > 0:
#           delete all Annotations for Annotation Type/Reference.
#       if mode = new and user != none :
#           delete all Annotations for Annotation Type/User.
#       if mode = new and reference = 0 and user = none:
#           delete all Annotations for Annotation Type.
#       if mode = append: do nothing special
#       if mode = preview:  set "DEBUG" to True
#
#    3. Load Evidence Codes and Terms into dictionaries for quicker lookup.
#           If not -O, then only load non-obsolete Terms.
#
#    4. For each line in the input file:
#
#       1.  Verify the Term Accession ID is valid for the given Annotation Type.
#           If the verification fails, report the error and skip the record.
#
#       2.  Verify the MGI ID is valid for the given Annotation Type.
#           If verification fails, report the error and skip the record.
#           If verification succeeds, store the MGI ID/Key pair in a dictionary
#           for future reference.
#
#       3.  Verify the J: is valid.
#           If verification fails, report the error and skip the record.
#           If verification succeeeds, store the Jnum/Key pair in a dictionary
#           for future reference.
#
#       4.  Verify the Evidence Code is valid for the given Annotation Type.
#           If verification fails, report the error and skip the record.
#
#       5.  Verify the Editor is provided (i.e. is not null).
#           If  verification fails, report the error and skip the record.
#
#       6.  If the Date is not given, use the current date.
#
#       7.  Determine if an Annotation record for the given Term, Object and
#           Qualifier already exists in VOC_Annot.
#               If it does exist, use its _Annot_key for the Evidence record.
#               If it does not exist, create a new Annotation record.
#
#       8.  Determine if an Evidence record for the given Annot key,
#           Evidence Code and Reference already exists in VOC_Evidence.
#               If it does exist, write a message to the error file.
#               If it does not exist, create a new Evidence record.
#
# History:
#
# lec   11/16/2023
#       wts2-1155/GOC taking over GOA mouse, GOA human, etc.
#
# sc	04/24/2017
#	- TR12556 - remove call to ALL_postMP
#
'''

import sys
import os
import re
import db
import accessionlib
import mgi_utils
import loadlib

libpath = os.environ['ANNOTLOAD'] + '/lib'
sys.path.insert(0, libpath)
import vocabloadlib

#db.setTrace(True)

# globals

# from configuration file
user = os.environ['PG_DBUSER']
passwordFileName = os.environ['PG_1LINE_PASSFILE']
mode = os.environ['ANNOTMODE']
inputFileName = os.environ['ANNOTINPUTFILE']
annotTypeName = os.environ['ANNOTTYPENAME']

annotProperty = 0
if 'ANNOTPROPERTY' in os.environ:
    annotProperty = os.environ['ANNOTPROPERTY']

delByReference = os.environ['DELETEREFERENCE']
delByUser = os.environ['DELETEUSER'] + '%'
loadObsolete = os.environ['ANNOTOBSOLETE']
loadType = 'Standard'	# the default load type of this annotload run

DEBUG = 0		# set DEBUG to false unless preview mode is selected

inputFile = ''		# file descriptor
diagFile = ''		# file descriptor
errorFile = ''		# file descriptor
annotFile = ''		# file descriptor
evidenceFile = ''	# file descriptor
propertyFile = ''	# file descriptor
noteFile = ''		# file descriptor

diagFileName = ''	# file name
errorFileName = ''	# file name
annotFileName = ''	# file name
evidenceFileName = ''	# file name
propertyFileName = ''	# file name
noteFileName = ''	# file name

delByReferenceKey = 0	# deletion reference key
annotTypeKey = 0	# VOC_AnnotType._AnnotType_key
annotKey = 0		# VOC_Annot._Annot_key
evidencePrimaryKey = 0	# VOC_Evidence._AnnotEvidence_key
propertyKey = 0		# VOC_Evidence_Property._EvidenceProperty_key
logicalDBKey = 1	# ACC_Accession._LogicalDB_key (default is "MGI", 1)
noteKey = 0		# MGI_Note._Note_key
mgiNoteObjectKey = 25	# MGI_Note._MGIType_key
mgiNoteTypeKey = 1008	# MGI_Note._NoteType_key

termDict = {}		# dictionary of terms for quick lookup
objectDict = {}		# dictionary of objects for quick lookup
referenceDict = {}	# dictionary of references for quick lookup
annotDict = {}		# dictionary of annotation records for quick lookup
evidenceDict = {}	# dictionary of evidence records for quick lookup
pTermDict = {}		# dictionary of propery terms for quick lookup

loaddate = loadlib.loaddate

# execute this SQL after the bcp files
execSQL = ''

# true (1) if this is the mcvload
isMCV = 0

# true (1) if this is a MP load
isMP = 0

# true (1) if this is a goload
isGO = 0

# true (1) if this is a disease/marker, disease/allele rollup load
isDiseaseMarker = 0
isDiseaseAllele = 0

# true (1) if this is an MP/marker, MP/allele rollup load
isMPMarker = 0
isMPAllele = 0

# true (1) if this is an OMIM/HPO load
isOMIMHPO = 0

# true (1) if no bcp files to load
skipBCP = 1

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
    global annotFile, annotFileName, evidenceFile, evidenceFileName
    global propertyFile, propertyFileName
    global noteFile, noteFileName
    global annotTypeKey, annotKey, annotTypeName, evidencePrimaryKey
    global noteKey, propertyKey
    global isMCV, isMP, isGO, isDiseaseMarker, isDiseaseAllele, isMPMarker, isMPAllele, isOMIMHPO
    global loadType

    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)

    head, tail = os.path.split(inputFileName) 

    diagFileName = tail + '.diagnostics'
    errorFileName = tail + '.error'
    annotFileName = tail + '.VOC_Annot.bcp'
    evidenceFileName = tail + '.VOC_Evidence.bcp'
    propertyFileName = tail + '.VOC_Evidence_Property.bcp'
    noteFileName = tail + '.MGI_Note.bcp'

    # determine load-type

    if len(sys.argv) == 2:

        loadType = sys.argv[1]
        print('LOAD TYPE: %s' % sys.argv[1])

        if loadType == 'mcv':
            isMCV = 1

        elif loadType == 'mp':
            isMP = 1

        elif loadType == 'go':
            isGO = 1

        elif loadType == 'diseaseMarker':
            isDiseaseMarker = 1
        
        elif loadType == 'diseaseAllele':
            isDiseaseAllele = 1

        elif loadType == 'mpMarker':
            isMPMarker = 1

        elif loadType == 'mpAllele':
            isMPAllele = 1

        elif loadType == 'omimhpo':
            isOMIMHPO = 1
        
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
        propertyFile = open(propertyFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % propertyFileName)
            
    try:
        noteFile = open(noteFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % noteFileName)
            
    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    # Set Log File Descriptor
    #db.set_sqlLogFD(diagFile)

    diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.write('Server: %s\n' % (db.get_sqlServer()))
    diagFile.write('Database: %s\n' % (db.get_sqlDatabase()))
    diagFile.write('Annotation Type Name: %s\n' % (annotTypeName))
    diagFile.write('Annotation File: %s\n' % (inputFileName))
    diagFile.write('Deletion Reference: %s\n\n' % (delByReference))
    diagFile.write('Deletion User: %s\n\n' % (delByUser))
    diagFile.write('Load Type: %s\n\n' % (loadType))

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

    results = db.sql(''' select _AnnotType_key from VOC_AnnotType where name = '%s' ''' % (annotTypeName), 'auto')

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

    global DEBUG, delByReferenceKey

    if delByReference != "J:0":

        delByReferenceKey = verifyReference(delByReference, 0, '')

        if delByReferenceKey is None:
            exit(1, 'Invalid Reference: %s\n' % (delByReference))

    if mode == 'new' or mode == 'delete':

        # verify deletion reference

        #
        # create temp table toDelete
        #
        if delByReference != "J:0":
            db.sql('''
                create temp table toDelete as 
                select e._Annot_key, e._AnnotEvidence_key 
                from VOC_Annot a, VOC_Evidence e
                where e._Refs_key = %s
                and e._Annot_key = a._Annot_key 
                and a._AnnotType_key = %s
                ''' % (delByReferenceKey, annotTypeKey), None, execute = not DEBUG)

        elif delByUser != "none%":
            # updated this to include the wildcard at the end of the comparison
            # above, as it is always added to the DELETEUSER environment
            # variable.  (so this code was always getting executed, even when
            # the config file had the DELETEUSER set to 'none') - jsb, 11/3/14

            db.sql('''
                create temp table toDelete as 
                select e._Annot_key, e._AnnotEvidence_key
                from VOC_Annot a, VOC_Evidence e, MGI_User u 
                where e._CreatedBy_key = u._User_key 
                and u.login like '%s'
                and e._Annot_key = a._Annot_key
                and a._AnnotType_key = %s
                ''' % (delByUser, annotTypeKey), None, execute = not DEBUG)

        else:
            db.sql('''
                create temp table toDelete as 
                select e._Annot_key, e._AnnotEvidence_key 
                from VOC_Annot a, VOC_Evidence e
                where a._AnnotType_key = %s
                and a._Annot_key = e._Annot_key 
                ''' % (annotTypeKey), None, execute = not DEBUG)

        db.sql('create index td_idx1 on toDelete(_Annot_key)', None)
        db.sql('create index td_idx2 on toDelete(_AnnotEvidence_key)', None)
            
        db.sql('''
            delete from MGI_Note n
            using toDelete d, VOC_Evidence_Property p
            where d._AnnotEvidence_key = p._AnnotEvidence_key
            and p._EvidenceProperty_key = n._Object_key 
            and n._MGIType_key = 41
            ''', None, execute = not DEBUG)

        db.sql('''
            delete from VOC_Evidence_Property p
            using toDelete d
            where d._AnnotEvidence_key = p._AnnotEvidence_key
            ''', None, execute = not DEBUG)

        db.sql('''
            delete from VOC_Evidence e
            using toDelete d
            where d._AnnotEvidence_key = e._AnnotEvidence_key
            ''', None, execute = not DEBUG)

        db.sql('''
            delete from VOC_Annot a
            using toDelete d
            where d._Annot_key = a._Annot_key
            and not exists (select 1 from VOC_Evidence e where d._Annot_key = e._Annot_key)
            ''', None, execute = not DEBUG)

        # remove the Used-FC references
        if isMP:
            db.sql('''
                delete from MGI_Reference_Assoc
                where _MGIType_key = 11
                and _Refs_key = %s
                and _RefAssocType_key = 1017
                ''' % (delByReferenceKey), None, execute = not DEBUG)

        db.sql('''drop table todelete''', None)
        db.commit()
            
    elif mode == 'append':
        pass

    elif mode == 'preview':
        DEBUG = 1
    else:
        exit(1, 'Invalid Processing Mode:  %s\n' % (mode))

    if mode == 'delete':
        exit(0)

def verifyTerm(termID, lineNum, line):
    '''
    # requires:
    #	termID - the Accession ID of the term
    #	lineNum - the line number of the record from the input file
    #   line = the line itself
    #
    # effects:
    #	verifies that:
    #		the Term exists 
    #		is of the appropriate type for the Annotation Type 
    #		    of the load 
    #	writes to the error file if the Term is invalid
    #
    # returns:
    #	0 if the Term is invalid
    #	Term Key if the Term is valid
    #
    '''

    termKey = 0

    if termID in termDict:
        termKey = termDict[termID]
    else:
        if loadObsolete == '0':
            errorFile.write('Invalid or Obsolete Term (%d): %s\n%s\n' % (lineNum, termID, line))
        else:
            errorFile.write('Invalid Term (%d): %s\n%s\n' % (lineNum, termID, line))
        termKey = 0

    return(termKey)

def verifyReference(referenceID, lineNum, line):
    '''
    # requires:
    #	referenceID - the Accession ID of the reference
    #	lineNum - the line number of the record from the input file
    #   line = the line itself
    #
    # effects:
    #	verifies that:
    #		the Reference  exists 
    #	writes to the error file if the Reference is invalid
    #
    # returns:
    #	0 if the Term is invalid
    #	Term Key if the Term is valid
    #
    '''

    referenceKey = 0

    if referenceID in referenceDict:
        referenceKey = referenceDict[referenceID]
    else:
        errorFile.write('Invalid Reference (%d): %s\n%s\n' % (lineNum, referenceID, line))
        referenceKey = 0

    return(referenceKey)

def verifyObject(objectID, logicalDBKey, lineNum, line):
    '''
    # requires:
    #	objectID - the ID of the Object
    #	logicalDBKey - the logical DB of the objectID
    #	lineNum - the line number of the record from the input file
    #   line = the line itself
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

    if objectID in objectDict:
            objectKey = objectDict[objectID]
    else:
        results = db.sql('''
            select a._Object_key 
            from ACC_Accession a, VOC_AnnotType t
            where lower(a.accID) = lower('%s')
            and a._LogicalDB_key = %s
            and a._MGIType_key = t._MGIType_key
            and t._AnnotType_key = %s
            ''' % (objectID, logicalDBKey, annotTypeKey), 'auto')

        if len(results) == 1:
            objectKey = results[0]['_Object_key']
            objectDict[objectID] = objectKey
        else:
            errorFile.write('Invalid Object (%d): %s\n%s\n' % (lineNum, objectID, line))
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

    global annotKey, evidencePrimaryKey, noteKey, propertyKey

    results = db.sql(''' select nextval('voc_annot_seq') as maxKey ''', 'auto')
    annotKey = results[0]['maxKey']

    results = db.sql(''' select nextval('voc_evidence_seq') as maxKey ''', 'auto')
    evidencePrimaryKey = results[0]['maxKey']

    results = db.sql(''' select nextval('mgi_note_seq') as maxKey ''', 'auto')
    noteKey = results[0]['maxKey']

    results = db.sql(''' select nextval('voc_evidence_property_seq') as maxKey ''', 'auto')
    propertyKey = results[0]['maxKey']

def loadReferenceDictionary():
    '''
    # requires:
    #
    # effects:
    #	loads global dictionaries (see global below) for quicker lookup
    #
    # returns:
    #	nothing
    '''

    global referenceDict

    #   
    # referenceDict
    #   
    results = db.sql('''
       select _Object_key, accID
       from ACC_Accession
       where _MGIType_key = 1 
       and _LogicalDB_key = 1 
       and prefixPart = 'J:'
       ''', 'auto')
    for r in results:
        referenceDict[r['accID']] = r['_Object_key']

def loadDictionaries():
    '''
    # requires:
    #
    # effects:
    #	loads global dictionaries (see global below) for quicker lookup
    #
    # returns:
    #	nothing
    '''

    global termDict, annotDict, evidenceDict, pTermDict

    # cache annotation type vocabulary

    cmd = '''
        select t._Object_key, t.accID 
        from VOC_Term_Acc_View t, VOC_Term tm, VOC_AnnotType a
        where t._Object_key = tm._Term_key 
        and tm._Vocab_key = a._Vocab_key
        and a._AnnotType_key = %s''' % (annotTypeKey)

    # if loadObsolete is false, then only load non-obsoleted terms...

    if loadObsolete == '0':
        # only load non-obsoleted terms
        cmd = cmd + ' and tm.isObsolete = 0'

    results = db.sql(cmd, 'auto')

    for r in results:
        termDict[r['accID']] = r['_Object_key']

    # cache property vocabulary(s)

    cmd = '''
        select _Term_key, term
        from VOC_Term
        where _Vocab_key in (%s)
        ''' % (annotProperty)
    results = db.sql(cmd, 'auto')

    for r in results:
        pTermDict[r['term']] = r['_Term_key']

    # cache annotation keys for this type of annotation

    results = db.sql('''
        select _Annot_key, _Object_key, _Term_key, _Qualifier_key 
        from VOC_Annot
        where _AnnotType_key = %s
        ''' % (annotTypeKey), 'auto')
    for r in results:
        key = '%s:%s:%s:%s' % (annotTypeKey, r['_Object_key'], r['_Term_key'], r['_Qualifier_key'])
        value = r['_Annot_key']
        annotDict[key] = value

    # cache evidence keys for this type of annotation

    cmd = '''
        select e._Annot_key, e._EvidenceTerm_key, e._Refs_key 
        from VOC_Evidence e, VOC_Annot a
        where a._AnnotType_key = %s
        and a._Annot_key = e._Annot_key''' % (annotTypeKey)

    results = db.sql(cmd, 'auto')
    for r in results:
        key = '%s:%s:%s' % (r['_Annot_key'], r['_EvidenceTerm_key'], r['_Refs_key'])
        value = r['_Annot_key']
        evidenceDict[key] = value

def loadObjectDict():
    global objectDict

    # cache object keys

    print('ldb: %s, annotType: %s' % (logicalDBKey, annotTypeKey))

    results = db.sql('''
        select a.accID, a._Object_key
        from ACC_Accession a, VOC_AnnotType t
        where a._LogicalDB_key = %s
        and a.preferred = 1
        and a._MGIType_key = t._MGIType_key
        and t._AnnotType_key = %s''' % (logicalDBKey, annotTypeKey), 'auto')
    for r in results:
        key = r['accID']
        value = r['_Object_key']
        objectDict[key] = value

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
    global execSQL

    # if an annotation already exists for the same 
    # AnnotType/Object/Term/Qualifier, use the same annotation key

    aKey = '%s:%s:%s:%s' % (annotTypeKey, objectKey, termKey, qualifierKey)

    # annotation may exist in our dictionary already...

    if aKey in annotDict:
        useAnnotKey = annotDict[aKey]
    else:
        useAnnotKey = annotKey
        annotDict[aKey] = useAnnotKey
        annotKey = annotKey + 1

        # create the new VOC_Annot record

        annotFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                % (useAnnotKey, annotTypeKey, objectKey, termKey, qualifierKey, entryDate, entryDate))

    return(useAnnotKey)

def createEvidenceRecord(newAnnotKey, evidenceKey, referenceKey, \
        inferredFrom, editorKey, notes, properties, entryDate, line, lineNum):
    '''
    # requires:
    #	newAnnotKey - primary key of the Annotation object
    #	evidenceKey - primary key of the Evidence Code
    #	referenceKey - primary key of the Reference
    #	inferredFrom - inferred from value
    #	editorKey - primary key of the Editor
    #	notes - notes
    #	properties - properties
    #	entryDate - creation and modification date of Annotation
    #	line - full line from input file
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

    global evidencePrimaryKey, evidenceDict, noteKey, propertyKey

    #
    # evidenceDict is used to check for duplicates:
    #	1) existing annotation (in the database) : see loadDictionaries/evidenceDict)
    #

    if isMP or isOMIMHPO:
            eKey = '%s:%s:%s:%s' % (newAnnotKey, evidenceKey, referenceKey, properties)

    elif isMPMarker or isMPAllele:
            eKey = '%s:%s:%s:%s:%s' % (newAnnotKey, evidenceKey, referenceKey, properties, notes)

    elif isGO or isDiseaseMarker or isDiseaseAllele:
            eKey = '%s:%s:%s:%s:%s' % (newAnnotKey, evidenceKey, referenceKey, properties, inferredFrom)

    # default duplication check
    else:
            eKey = '%s:%s:%s' % (newAnnotKey, evidenceKey, referenceKey)

    # evidence record may exist in our dictionary already
    # if so, it's a duplicate; let's report it
    if eKey in evidenceDict:
            errorFile.write('Duplicate evidence (%d): \n%s\n' % (lineNum, line))
            return

    else:
        # not a duplicate
        evidenceDict[eKey] = eKey

    evidenceFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            % (evidencePrimaryKey, newAnnotKey, evidenceKey, referenceKey, \
                inferredFrom, editorKey, editorKey, entryDate, entryDate))

    # storing data in MGI_Note
    if len(notes) > 0:

        noteFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            % (noteKey, evidencePrimaryKey, mgiNoteObjectKey, mgiNoteTypeKey, notes, \
               editorKey, editorKey, entryDate, entryDate))

        noteKey = noteKey + 1

    #
    # storing data in VOC_Evidence_Property
    #
    # for each property in the property field 11
    #     extract the stanza
    #       extract the property/value pair
    #         extract the term and the value
    #           if the term exists in the pTermDict lookup then
    #             create a row in the property table
    #

    if len(properties) > 0:

        stanza = 1
        allStanzas = str.split(properties, '&===&')

        for s in allStanzas:

            seqnum = 1
            allProps = str.split(s, '&==&')

            for p in allProps:

                #printproperties
                pTerm, pValue = str.split(p,'&=&')

                if pTerm in pTermDict:
                    propertyFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                    % (propertyKey, evidencePrimaryKey, pTermDict[pTerm], stanza, seqnum, pValue, \
                       editorKey, editorKey, entryDate, entryDate))

                    seqnum = seqnum + 1
                    propertyKey = propertyKey + 1
                else:
                    errorFile.write('Invalid Property (%s):  %s\n%s\n' % (lineNum, pTerm, line))

            stanza = stanza + 1

    evidencePrimaryKey = evidencePrimaryKey + 1

def processMcvFile():
    '''
    # requires:
    #
    # effects:
    #       Reads Marker Category input file
    #       Verifies and Processes each line in the input file
    # 	Deletes all annotations to all markers in the input file
    #	Creates new annotations  in the input file, unless a delete is
    #	indicated (all fields blank except mgi ID)
    #
    # returns:
    #       nothing
    #
    '''
    global logicalDBKey, annotDict, skipBCP

    skipBCP = 0

    # running list of markers in the input so we don't delete more than once
    mkrKeyList = []

    #  we first delete all annotations by marker, 
    # then recreate based on the input
    deleteSQL = '''delete from VOC_Annot where _AnnotType_key = %s and _Object_key = %s'''
    lineNum = 0

    # For each line in the input file
    for line in inputFile.readlines():
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = str.split(line[:-1], '\t')

        try:
            termID = tokens[0]
            mgiID = tokens[1]
            jnum = tokens[2]
            evidence = tokens[3]
            inferredFrom = str.strip(tokens[4])
            qualifier = str.strip(tokens[5])
            editor = str.strip(tokens[6])
            entryDate = str.strip(tokens[7])
            notes = str.strip(tokens[8])
            properties = ''

            if len(tokens) > 9:
                # field 10 reserved for optional ldb the default is "1" (MGI)
                col10 = accessionlib.get_LogicalDB_key(tokens[9])
                if col10 != None:
                    logicalDBKey = col10
            if len(tokens) > 10:
                # field 11 reserved for optional properties
                properties = str.strip(tokens[10])

        except:
            exit(1, 'Invalid Line (%d): \n%s\n' % (lineNum, line))

        # if we just have an MGI ID delete all annotations and continue
        if termID == '' and jnum == '' and evidence == '' and \
           qualifier == '' and editor == '':

            markerKey = verifyObject(mgiID, logicalDBKey, lineNum, line)

            if markerKey == 0: # not valid
                # skip this record, verifyObject logs to discrepancy file
                print('skipping record, invalid MGI ID for delete only mode')
                continue
            else:
                # delete existing annotations for this marker/annotation type
                # and go on to next record
                print('delete only mode, deleting all annotations for %s' % markerKey)
                print(deleteSQL % (annotTypeKey, markerKey))
                db.sql(deleteSQL % (annotTypeKey, markerKey), None)
                continue

        # if we get here, continue verifying
        termKey = verifyTerm(termID, lineNum, line)
        markerKey = verifyObject(mgiID, logicalDBKey, lineNum, line)
        referenceKey = verifyReference(jnum, lineNum, line)
        evidenceKey = vocabloadlib.verifyEvidence(evidence, annotTypeKey, lineNum, errorFile)
        qualifierKey = vocabloadlib.verifyQualifier(qualifier, annotTypeKey, 0, lineNum, errorFile)
        editorKey = loadlib.verifyUser(editor, lineNum, errorFile)
        
        # if any verification failed, this is an error
        if termKey == 0 or markerKey == 0 or \
           referenceKey == 0 or \
           evidenceKey == 0 or \
           qualifierKey == 0 or \
           editorKey == 0:

            continue

        # If the entry date is not given, use the current date
        if len(entryDate) == 0:
            entryDate = loaddate

        # if we get here, there are no errors so process the annotation

        # first delete if we haven't already seen this marker in the input
        if not markerKey in mkrKeyList:
            db.sql(deleteSQL % (annotTypeKey, markerKey), None)
            mkrKeyList.append(markerKey)

        # then create annotations
        # first delete from annotDict, because we've deleted from the database
        aKey = '%s:%s:%s:%s' % (annotTypeKey, markerKey, termKey, qualifierKey)
        if aKey in annotDict:
            del annotDict[aKey]

        # now create an annotation record
        newAnnotKey = createAnnotationRecord( \
            markerKey, termKey, qualifierKey, entryDate)
        
        # now create evidence record
        createEvidenceRecord(newAnnotKey, \
            evidenceKey, \
            referenceKey, \
            inferredFrom, \
            editorKey, \
            notes, \
            properties, \
            entryDate, \
            line, \
            lineNum)

    db.commit()

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
    global skipBCP

    lineNum = 0

    # For each line in the input file
    objectDictLoaded = 0

    for line in inputFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = str.split(line[:-1], '\t')

        try:
            termID = tokens[0]
            objectID = tokens[1]
            jnum = tokens[2]
            evidence = tokens[3]
            inferredFrom = str.strip(tokens[4])
            qualifier = str.strip(tokens[5])
            editor = str.strip(tokens[6])
            entryDate = str.strip(tokens[7])
            notes = str.strip(tokens[8])
            properties = ''

            if len(tokens) > 9:
                # field 10 reserved for optional ldb
                # the default is "1" (MGI)
                col10 = accessionlib.get_LogicalDB_key(str.strip(tokens[9]))

                if col10 != None:
                    logicalDBKey = col10

                if len(tokens) > 10:
                    # field 11 reserved for optional properties
                    properties = str.strip(tokens[10])

        except:
            exit(1, 'Invalid Line (%d): \n%s\n' % (lineNum, line))

        # for files that specify logicalDB - we won't know that value until we get the first term
        if objectDictLoaded == 0:
            loadObjectDict()
            objectDictLoaded = 1

        termKey = verifyTerm(termID, lineNum, line)
        objectKey = verifyObject(objectID, logicalDBKey, lineNum, line)
        referenceKey = verifyReference(jnum, lineNum, line)
        evidenceKey = vocabloadlib.verifyEvidence(evidence, annotTypeKey, lineNum, errorFile)
        qualifierKey = vocabloadlib.verifyQualifier(qualifier, annotTypeKey, 0, lineNum, errorFile)
        editorKey = loadlib.verifyUser(editor, lineNum, errorFile)

        if termKey == 0 or \
            objectKey == 0 or \
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
            properties, \
            entryDate, \
            line, \
            lineNum)

        # don't skip the bcp file loading...data exists that needs to be loaded
        skipBCP = 0

    # end of "for line in inputFile.readlines():"
    db.commit()

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

    global execSQL

    annotFile.close()
    evidenceFile.close()
    noteFile.close()
    propertyFile.close()

    if DEBUG:
        return

    if skipBCP:
        return 0

    db.commit()

    bcpCommand = os.environ['PG_DBUTILS'] + '/bin/bcpin.csh'
    currentDir = os.getcwd()

    bcpAnnotCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
        (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(),'VOC_Annot',
        currentDir, annotFileName)
    diagFile.write('%s\n' % bcpAnnotCmd)

    bcpEvidenceCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
        (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(),'VOC_Evidence',
        currentDir, evidenceFileName)
    diagFile.write('%s\n' % bcpEvidenceCmd)

    bcpNoteCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
        (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(),'MGI_Note',
        currentDir, noteFileName)
    diagFile.write('%s\n' % bcpNoteCmd)

    bcpPropertyCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
        (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(),'VOC_Evidence_Property',
        currentDir, propertyFileName)
    diagFile.write('%s\n' % bcpPropertyCmd)

    print('BCPing files')
    os.system(bcpAnnotCmd)
    os.system(bcpEvidenceCmd)
    os.system(bcpNoteCmd)
    os.system(bcpPropertyCmd)
    print('BCP done')

    # update voc_annot_seq auto-sequence
    execSQL = '''select setval('voc_annot_seq', (select max(_Annot_key) from VOC_Annot))'''
    print(execSQL)
    db.sql(execSQL, None)
    db.commit()

    # update voc_evidence_seq auto-sequence
    execSQL = '''select setval('voc_evidence_seq', (select max(_AnnotEvidence_key) from VOC_Evidence))'''
    print(execSQL)
    db.sql(execSQL, None)
    db.commit()

    # update voc_evidence_property_seq auto-sequence
    execSQL = '''select setval('voc_evidence_property_seq', (select max(_EvidenceProperty_key) from VOC_Evidence_Property))'''
    print(execSQL)
    db.sql(execSQL, None)
    db.commit()

    # update mgi_note_seq auto-sequence
    execSQL = '''select setval('mgi_note_seq', (select max(_Note_key) from MGI_Note))'''
    print(execSQL)
    db.sql(''' select setval('mgi_note_seq', (select max(_Note_key) from MGI_Note)) ''', None)
    db.commit()

#
# Main
#

print('\nannotload.py - main() started')

#print('\nannotload.py - init')
init()

#print('\nannotload.py - verifyAnnotType')
verifyAnnotType()

#print('\nannotload.py - loadReferenceDictionary')
loadReferenceDictionary()

#print('\nannotload.py - verifyMode')
verifyMode()

#print('\nannotload.py - setPrimaryKeys')
setPrimaryKeys()

#print('\nannotload.py - loadDictionary')
loadDictionaries()

#print('\nannotload.py - process')
if isMCV:
    processMcvFile()
else:
    processFile()

bcpFiles()

print('\nannotload.py - main() finished')
