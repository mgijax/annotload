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
setenv ANNOTATIONFILE		${INPUTFILE}.annotload
setenv DELETEREFERENCE		"J:0"

setenv DBUSER			mgd_dbo
setenv DBPASSWORDFILE		${DBUTILITIESPATH}/.mgd_dbo_password

cd `dirname $0`

# create the Annotation File
gomarker.py -S${DBSERVER} -D${DBNAME} -I${INPUTFILE}

# load the Annotation File
annotload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${MODE} -I${ANNOTATIONFILE} -A\"${ANNOTATIONTYPENAME}\" -R${DELETEREFERENCE}
