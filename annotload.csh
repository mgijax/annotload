#!/bin/csh -f -x

#
# Annotation Load Wrapper
#
# Usage:
# 	annotload.csh config file
#

setenv CONFIGFILE $1

cd `dirname 0` && source ./${CONFIGFILE}

cd ${ANNOTDATADIR}
${ANNOTLOAD} -S${MGD_DBSERVER} -D${MGD_DBNAME} -U${MGI_DBUSER} -P${MGI_DBPASSWORDFILE} -M${ANNOTMODE} -I${ANNOTFILE} -A"${ANNOTTYPENAME}" -R${DELETEREFERENCE}
