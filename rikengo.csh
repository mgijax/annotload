#!/bin/csh -f -x

#
# RIKEN GO Annotation Load Wrapper
#
# Usage:
# 	rikengo.csh DBSERVER DBNAME inputfile mode
#
# Purpose:
#	executes rikengo.py (to generate a list of new genes)
#	adds the new genes to Nomen and broadcasts them to MGD
#	executes rikengo.py again (to generate the annotation loader input file)
#	executes annotload.py
#
#	To be used only in a development database (for now).
#

setenv DBSERVER		$1
setenv DBNAME		$2
setenv INPUTFILE	$3
setenv MODE		$4

setenv DBUTILITIESPATH		/usr/local/mgi/dbutils/mgidbutilities
setenv ANNOTATIONTYPENAME	"GO/Marker"
setenv ANNOTATIONFILE		`basename ${INPUTFILE}`.annotload
setenv DELETEREFERENCE		"J:23000"

setenv DBUSER			mgd_dbo
setenv DBPASSWORDFILE		${DBUTILITIESPATH}/.mgd_dbo_password

setenv SCHEMADIR		/home/lec/db/mgddbschema
setenv PREPUBLOAD		/home/lec/loads/rikenprepubload
setenv POSTPUBLOAD		/home/lec/loads/rikenpostpubload

echo 'RIKEN GO Annotation Load'
date

set annotdir = `dirname $0`

# create the new genes file
${annotdir}/rikengo.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -I${INPUTFILE}

# process the new genes file (add new genes to nomen and broadcast them to mgd)
cp ${INPUTFILE}.annotload.newgenes ${PREPUBLOAD}/data
rm -rf ${PREPUBLOAD}/data/new_genes.txt
ln -s ${PREPUBLOAD}/data/${INPUTFILE}.annotload.newgenes ${PREPUBLOAD}/data/new_genes.txt
${PREPUBLOAD}/addToNomen.sh
${POSTPUBLOAD}/addGenes.sh

# truncate the annotation tables
${SCHEMADIR}/table/VOC_Annot_truncate.object
${SCHEMADIR}/table/VOC_Evidence_truncate.object

# drop indexes
${SCHEMADIR}/index/VOC_Annot_drop.object
${SCHEMADIR}/index/VOC_Evidence_drop.object

# create the Annotation File
${annotdir}/rikengo.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -I${INPUTFILE}

# load the Annotation File; load annotations to obsolete terms
${annotdir}/annotload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${MODE} -I${ANNOTATIONFILE} -A\"${ANNOTATIONTYPENAME}\" -R${DELETEREFERENCE} -O

# create indexes
${SCHEMADIR}/index/VOC_Annot_create.object
${SCHEMADIR}/index/VOC_Evidence_create.object

date

