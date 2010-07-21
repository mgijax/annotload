#!/bin/csh -f -x

#
# Annotation Load Wrapper
#
# Usage:
# 	annotload.csh config file
#

setenv CONFIGFILE $1
setenv LOADTYPE $2

source ${CONFIGFILE}

echo 'running annotload.csh'
rm -rf ${ANNOTLOG}
touch ${ANNOTLOG}

if ( ! -r ${ANNOTINPUTFILE} ) then
    echo "Cannot read input file: ${ANNOTINPUTFILE}"
    exit 1
endif

date >> ${ANNOTLOG}
${ANNOTLOAD}/annotload.py ${LOADTYPE} | tee -a ${ANNOTLOG}
date >> ${ANNOTLOG}
