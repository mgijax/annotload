#!/bin/csh -f -x

#
# PhenoSlim/Genotype Annotation Load Wrapper
#
# Usage:
# 	phenoslimgenotype.csh DBSERVER DBNAME inputfile editor mode
#
# Purpose:
#	executes phenoslimgenotype.py
#	executes annotload.py if the annotation file exists
#

setenv DBSERVER		$1
setenv DBNAME		$2
setenv INPUTFILE 	$3
setenv EDITOR 		$4
setenv MODE 		$5

setenv DBUTILITIESPATH		/usr/local/mgi/dbutils/mgidbutilities
setenv ANNOTATIONTYPENAME	"PhenoSlim/Genotype"
setenv ANNOTATIONFILE		${INPUTFILE}.annotload
setenv DELETEREFERENCE		"J:0"

setenv DBUSER			mgd_dbo
setenv DBPASSWORDFILE		${DBUTILITIESPATH}/.mgd_dbo_password

echo 'PhenoSlim/Genotype Annotation Load'
date

cd `dirname $0`

# create the Annotation File
phenoslimgenotype.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -I${INPUTFILE} -E{$EDITOR}

# load the Annotation File
annotload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${MODE} -I${ANNOTATIONFILE} -A\"${ANNOTATIONTYPENAME}\" -R${DELETEREFERENCE}

date

