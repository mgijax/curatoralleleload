#!/bin/sh
#
#  alleleQC.sh
###########################################################################
#
#  Purpose:
#
#      This script is a wrapper around the process 
#	that does QC checks for the Curator Allele load
#
#  Usage:
#
#      alleleQC.sh  filename  
#
#      where
#          filename = full path to the input file
#
#  Env Vars:
#
#      See the configuration file
#
#  Inputs:
#	Allele file
#
#  Outputs:
#
#      - QC report for the input file 	
#      - Log file (${QC_LOGFILE})
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  Sanity error (prompt to view sanity report)
#      2:  Unexpected error occured running alleleQC.py (prompt to view log)
#      3:  QC errors (prompt to view qc report)
#
#  Assumes:  Nothing
#
#  Implementation:
#
#      This script will perform following steps:
#
#      1) Validate the arguments to the script
#      2) Validate & source the configuration files to establish the environment
#      3) Verify that the input file exists
#      4) Update path to QC reports if this is not a 'live' run 
#	     i.e. curators running the scripts 
#      5) Initialize the log file
#      6) creates table in tempdb for the input file
#      7) Call alleleQC.py to generate the QC report
#      8) drops the tempdb table
#
#
#  Notes:  None
#
###########################################################################
#
#  Modification History:
#
#  Date        SE   Change Description
#  ----------  ---  -------------------------------------------------------
#
#  03/11/2014  sc  Initial development
#
###########################################################################
CURRENTDIR=`pwd`
BINDIR=`dirname $0`

CONFIG=`cd ${BINDIR}/..; pwd`/curatoralleleload.config
USAGE='Usage: alleleQC.sh  filename'

# set LIVE_RUN  to QC check only as the default
LIVE_RUN=0; export LIVE_RUN

#
# Make sure an input file was passed to the script. If the optional "live"
# argument is given, that means that the output files are located in the
# /data/loads/... directory, not in the current directory.
#
if [ $# -eq 1 ]
then
    INPUT_FILE=$1
elif [ $# -eq 2 -a "$2" = "live" ]
then
    INPUT_FILE=$1
    LIVE_RUN=1
else
    echo ${USAGE}; exit 1
fi

#echo "INPUT_FILE: ${INPUT_FILE}"

#
# Make sure the configuration file exists and source it.
#
if [ -f ${CONFIG} ]
then
    . ${CONFIG}
else
    echo "Missing configuration file: ${CONFIG}"
    exit 1
fi

#
# If the QC check is being run by a curator, the mgd_dbo password needs to
# be in a password file in their HOME directory because they won't have
# permission to read the password file in the pgdbutilities product.
#
if [ "${USER}" != "mgiadmin" ]
then
    PGPASSFILE=$HOME/.pgpass
fi

#
# If this is not a "live" run, the output, log and report files should reside
# in the current directory, so override the default settings.
#
if [ ${LIVE_RUN} -eq 0 ]
then
	QC_RPT=${CURRENTDIR}/`basename ${QC_RPT}`
	QC_LOGFILE=${CURRENTDIR}/`basename ${QC_LOGFILE}`

fi

#
# Initialize the log file.
#
LOG=${QC_LOGFILE}
rm -rf ${LOG}
touch ${LOG}

#
# Convert the input file into a QC-ready version that can be used to run
# the QC reports against.
#
dos2unix ${INPUT_FILE} ${INPUT_FILE} 2>/dev/null

#
# Create a temporary file and make sure it is removed when this script
# terminates.
#
TMP_FILE=/tmp/`basename $0`.$$
trap "rm -f ${TMP_FILE}" 0 1 2 15

#
# Make sure the input files exist (regular file or symbolic link).
#
if [ "`ls -L ${INPUT_FILE} 2>/dev/null`" = "" ]
then
    echo "" | tee -a ${LOG}
    echo "Input file does not exist: ${INPUT_FILE}" | tee -a ${LOG}
    echo "" | tee -a ${LOG}
    exit 1
fi

#
# Generate the QC reports.
#
echo "" >> ${LOG}
date >> ${LOG}
echo "Generate the QC reports" >> ${LOG}
{ ${PYTHON} ${CURATORALLELELOAD}/bin/alleleQC.py ${INPUT_FILE} 2>&1; echo $? > ${TMP_FILE}; } >> ${LOG}

if [ `cat ${TMP_FILE}` -eq 1 ]
then
    echo "An error occurred while generating the QC reports"
    echo "See log file (${LOG})"
    RC=1
elif [ `cat ${TMP_FILE}` -eq 2 ]
then
    if [ ${LIVE_RUN} -eq 0 ]
    then
	echo ""
	echo "Warn and Skipped Alleles in QC. See ${QC_RPT} " | tee -a ${LOG}
    fi
    RC=2
elif [ `cat ${TMP_FILE}` -eq 3 ]
then
    if [ ${LIVE_RUN} -eq 0 ]
    then
        echo ""
        echo "Skipped Alleles in QC. See ${QC_RPT} " | tee -a ${LOG}
    fi
    RC=3
elif [ `cat ${TMP_FILE}` -eq 4 ]
then
    if [ ${LIVE_RUN} -eq 0 ]
    then
        echo ""
        echo "Warn Alleles in QC. See ${QC_RPT} " | tee -a ${LOG}
    fi
    RC=4


else
    if [ ${LIVE_RUN} -eq 0 ]
    then
	echo "No QC errors detected"
    fi
    RC=0
fi

echo "" >> ${LOG}
date >> ${LOG}
echo "Finished running QC checks on the input file" >> ${LOG}

exit ${RC}
