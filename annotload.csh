#!/bin/csh -f -x

#
# Annotation Load Wrapper
#
# Usage:
# 	annotload.csh config file
#

setenv CONFIGFILE $1

source ${CONFIGFILE}

rm -rf ${ANNOTLOG}
touch ${ANNOTLOG}

date >> ${ANNOTLOG}

${ANNOTLOAD}/annotload.py | tee -a ${ANNOTLOG}

#
# delete/reload inferredfrom cache 
# for specific "deleted-by" user
#
# this call will be added as part of the new uniproload
# and will replace the independent calls to "inferredfrom"
# from individual loads (goaload, goratload, swissload)
#
#${MGICACHELOAD}/inferredfrom.py -S${MGD_DBSERVER} -D${MGD_DBNAME} -U${MGD_DBUSER} -P${MGD_DBPASSWORDFILE} -B"${DELETEUSER}" >> ${ANNOTLOG}

date >> ${ANNOTLOG}
