#!/bin/csh -f -x

#
# GO/Marker Annotation Load Wrapper
#
# Usage:
# 	gomarker.csh DBSERVER DBNAME inputfile mode
#
# Purpose:
#	executes gomarker.py
#	executes annotload.py
#

setenv DBSERVER		$1
setenv DBNAME		$2
setenv INPUTFILE	$3
setenv MODE		$4

setenv DBUTILITIESPATH		/usr/local/mgi/dbutils/mgidbutilities
setenv ANNOTATIONTYPENAME	"GO/Marker"
setenv ANNOTATIONFILE		`basename ${INPUTFILE}`.annotload
setenv DELETEREFERENCE		"J:0"

setenv DBUSER			mgd_dbo
setenv DBPASSWORDFILE		${DBUTILITIESPATH}/.mgd_dbo_password

echo 'GO/Marker Annotation Load'
date

cd `dirname $0`

# create the Annotation File
gomarker.py -S${DBSERVER} -D${DBNAME} -I${INPUTFILE}
# sort it by column 2 (MGI Marker ID)
sort -k 2,3,1,2 ${ANNOTATIONFILE} > ${ANNOTATIONFILE}.sorted
mv -f ${ANNOTATIONFILE}.sorted ${ANNOTATIONFILE}

# load the Annotation File
annotload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${MODE} -I${ANNOTATIONFILE} -A\"${ANNOTATIONTYPENAME}\" -R${DELETEREFERENCE}

date

