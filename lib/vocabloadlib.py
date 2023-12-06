#
# Program: vocabloadlib.py
# Purpose: Provide a set of library functions used by annotload
#

import sys
import os
import db

#globals

ecodeDict = {}         # evidence codes
qualifierDict = {}     # qualifiers

# Purpose:  verify Evidence Code
# Returns:  Evidence Code key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Evidence Code exists in the ecode dictionary
#	writes to the error file if the Evidence Code is invalid
# Throws:  nothing

def verifyEvidence(
    ecode, 	  # Evidence Code value (str.
    annotTypeKey, # Annotation Type key, (str.
    lineNum,	  # line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global ecodeDict

    ecodeKey = 0

    if len(ecodeDict) == 0:
        results = db.sql('''
                select e._Term_key, e.abbreviation
                from VOC_Term e, VOC_AnnotType t
                where e._Vocab_key = t._EvidenceVocab_key
                and t._AnnotType_key = %s
                ''' % (annotTypeKey), 'auto')

        for r in results:
            ecodeDict[r['abbreviation']] = r['_Term_key']

    if ecode in ecodeDict:
        ecodeKey = ecodeDict[ecode]
    else:
        if errorFile != None:
            errorFile.write('Invalid Evidence Code (%d): %s\n' % (lineNum, ecode))

    return ecodeKey

# Purpose:  verify Qualifier
# Returns:  Qualifier key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Qualifier exists in the qualifier dictionary
#	writes to the error file if the Qualifier is invalid
# Throws:  nothing

def verifyQualifier(
    qualifier, 	  # Qualifier value (str.
    annotTypeKey, # Annotation Type key, (str.
    byAbbrev, 	  # Compare by Abbreviation (1) or by Term (0)
    lineNum,	  # line number (integer)
    errorFile	  # error file (file descriptor)
    ):

    global qualifierDict

    qualifierKey = 0

    if len(qualifierDict) == 0:

        if byAbbrev:
            results = db.sql('''
                select e._Term_key, e.abbreviation
                from VOC_Term e, VOC_AnnotType t
                where e._Vocab_key = t._QualifierVocab_key 
                and t._AnnotType_key = %s
                ''' % (annotTypeKey), 'auto')

            for r in results:
                qualifierDict[r['abbreviation']] = r['_Term_key']

        else:
            results = db.sql('''
                select e._Term_key, e.term
                from VOC_Term e, VOC_AnnotType t
                where e._Vocab_key = t._QualifierVocab_key
                and t._AnnotType_key = %s
                ''' % (annotTypeKey), 'auto')

            for r in results:
                qualifierDict[r['term']] = r['_Term_key']

    if len(qualifier) == 0:
        qualifier = None

    if qualifier in qualifierDict:
        qualifierKey = qualifierDict[qualifier]
    else:
        if errorFile != None:
            errorFile.write('Invalid Qualifier (%d): %s\n' % (lineNum, qualifier))

    return qualifierKey
