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
date >> ${ANNOTLOG}
