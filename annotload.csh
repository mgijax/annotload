#!/bin/csh -f -x

#
# Annotation Load Wrapper
#
# Usage:
# 	annotload.csh DBSERVER DBNAME inputfile mode
#
# Purpose:
#	executes annotload.py
#

setenv SCHEMADIR	$1
setenv INPUTFILE	
setenv MODE		append
setenv ANNOTATIONTYPENAME	"Mammalian Phenotype/Genotype"
setenv ANNOTATIONFILE		`basename ${INPUTFILE}`.annotload
setenv DELETEREFERENCE		"J:0"

source ${SCHEMADIR}/Configuration

setenv ANNOTLOAD	/home/lec/loads/annotload-mpr

echo 'Annotation Load'
date

set annotdir = `dirname $0`

# load the Annotation File
${ANNOTLOAD}/annotload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${MODE} -I${ANNOTATIONFILE} -A"${ANNOTATIONTYPENAME}" -R${DELETEREFERENCE}

date

