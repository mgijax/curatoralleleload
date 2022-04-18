#!/bin/sh

#
# This script is a wrapper around the process that loads 
# Curated Alleles
#
#
#     curatoralleleload.sh 
#

cd `dirname $0`/..
CONFIG_LOAD=`pwd`/curatoralleleload.config

cd `dirname $0`
LOG=`pwd`/curatoralleleload.log
rm -rf ${LOG}

USAGE='Usage: curatoralleleload.sh'
SCHEMA='mgd'

#
#  Verify the argument(s) to the shell script.
#
if [ $# -ne 0 ]
then
    echo ${USAGE} | tee -a ${LOG}
    exit 1
fi

#
# verify & source the configuration file
#

if [ ! -r ${CONFIG_LOAD} ]
then
    echo "Cannot read configuration file: ${CONFIG_LOAD}"
    exit 1
fi

. ${CONFIG_LOAD}

#
# Just a verification of where we are at
#

echo "MGD_DBSERVER: ${MGD_DBSERVER}"
echo "MGD_DBNAME: ${MGD_DBNAME}"

#
#  Source the DLA library functions.
#

if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}" | tee -a ${LOG}
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined." | tee -a ${LOG}
    exit 1
fi

#
# verify input file exists and is readable
#

if [ ! -r ${INPUT_FILE_DEFAULT} ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "Cannot read from input file: ${INPUT_FILE_DEFAULT}"
fi

#
# createArchive including OUTPUTDIR, startLog, getConfigEnv
# sets "JOBKEY"
#

preload ${OUTPUTDIR}

#
# rm all files/dirs from OUTPUTDIR
#

cleanDir ${OUTPUTDIR}

# NOTE: keep this commented out until production release
#
# There should be a "lastrun" file in the input directory that was created
# the last time the load was run for this input file. If this file exists
# and is more recent than the input file, the load does not need to be run.
#
LASTRUN_FILE=${INPUTDIR}/lastrun
#if [ -f ${LASTRUN_FILE} ]
#then
#    if test ${LASTRUN_FILE} -nt ${INPUT_FILE_DEFAULT}
#    then
#
#        echo "Input file has not been updated - skipping load" | tee -a ${LOG_PROC}
#        # set STAT for shutdown
#        STAT=0
#        echo 'shutting down'
#        shutDown
#        exit 0
#    fi
#fi

echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo "Run QC checks"  | tee -a ${LOG_DIAG}
${CURATORALLELELOAD}/bin/alleleQC.sh ${INPUT_FILE_DEFAULT} live
STAT=$?
if [ ${STAT} -eq 1 ]
then
    checkStatus ${STAT} "An error occurred while generating the QC reports - See ${QC_LOGFILE}. alleleQC.sh"
    # run postload cleanup and email logs
    shutDown
fi

if [ ${STAT} -eq 2 ]
then
    echo "Warn and Skipped Alleles in QC. Warn alleles were loaded - See ${QC_RPT}. alleleQC.sh" | tee -a ${LOG_PROC} ${LOG_DIAG}


fi

if [ ${STAT} -eq 3 ]
then
    echo "Skipped Alleles in QC See ${QC_RPT}. alleleQC.sh" | tee -a ${LOG_PROC} ${LOG_DIAG}
    
fi

if [ ${STAT} -eq 4 ]
then
    echo "Warn Alleles in QC. These were loaded. See ${QC_RPT}. alleleQC.sh" | tee -a ${LOG_PROC} ${LOG_DIAG}

fi


#
# run the load
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo "Run curatoralleleload.py"  | tee -a ${LOG_DIAG}
${PYTHON} ${CURATORALLELELOAD}/bin/curatoralleleload.py  
STAT=$?
checkStatus ${STAT} "${CURATORALLELELOAD}/bin/curatoralleleload.py"

#
# Archive a copy of the input file, adding a timestamp suffix.
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo "Archive input file" >> ${LOG_DIAG}
TIMESTAMP=`date '+%Y%m%d.%H%M'`
# The published file
ARC_FILE=`basename ${INPUT_FILE_DEFAULT}`.${TIMESTAMP}
cp -p ${INPUT_FILE_DEFAULT} ${ARCHIVEDIR}/${ARC_FILE}

# the QC'd file
ARC_FILE=`basename ${INPUT_FILE_QC}`.${TIMESTAMP}
cp -p ${INPUT_FILE_QC} ${ARCHIVEDIR}/${ARC_FILE}

#
# Touch the "lastrun" file to note when the load was run.
#
if [ ${STAT} = 0 ]
then
    touch ${LASTRUN_FILE}
fi


# run postload cleanup and email logs

shutDown

