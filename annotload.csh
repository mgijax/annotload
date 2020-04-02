#!/bin/csh -f

#
# Annotation Load Wrapper
#
# Assumes: caller of this wrapper has sourced the
# master config file
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
${PYTHON} ${ANNOTLOAD}/annotload.py ${LOADTYPE} | tee -a ${ANNOTLOG}
set resultcode=$?
date >> ${ANNOTLOG}
exit $resultcode
