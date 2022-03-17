#
#  alleleQC.py
###########################################################################
#
#  Purpose:
#
#	This script will generate a QC report for a curator allele
#	    input file
#
#  Usage:
#
#      alleleQC.py  filename
#
#      where:
#          filename = path to the input file
#
#  Inputs:
#      - input file as parameter - see USAGE
#
#  Outputs:
#
#      - QC report (${QC_RPT})
#      - Warning report (${WARNING_RPT})
#      - temp table BCP file (${MGI_ID_BCP})
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#      2:  Fatal QC errors detected and written to report
#
#  Assumes:
#
#  Implementation:
#
#      This script will perform following steps:
#
#      1) Validate the arguments to the script.
#      2) Perform initialization steps.
#      3) Open the input/output files.
#      4) Generate the QC reports.
#      5) Close the input/output files.
#
#  History:
#
# 02/09/2022   sc   
#       - CREAM Project
#
###########################################################################

import sys
import os
import string
import db
import time

#
#  CONSTANTS
#
TAB = '\t'
CRT = '\n'

USAGE = 'Usage: alleleQC.py  inputFile'

#
#  GLOBALS
#

# Report file names
qcRptFile = os.environ['QC_RPT']
warnRptFile = os.environ['WARNING_RPT']

# 1 if any QC errors in the input file
hasFatalErrors = 0
hasWarnErrors = 0

# for bcp
bcpin = '%s/bin/bcpin.csh' % os.environ['PG_DBUTILS']
server = os.environ['MGD_DBSERVER']
database = os.environ['MGD_DBNAME']

# Lookups
geneIdLookup = []
userLookup = []
statusLookup = []
typeLookup = []
inheritModeLookup = []
transmissionLookup = []
collectionLookup = []
referenceLookup = []
pclLookup = []
strainLookup = []
mclLookup = []
subtypeLookup = []
mutationLookup = []

# report lists
dupeLineList = []
missingColumnList = []
reqColumnList = []
badGeneIdList = []
badUserList = []
badTypeList = []
badStatusList = []
badInheritModeList = []
imOSNnoGenNoteList = [] # inheritance mode is 'Other (see notes)' and there is not general note
badTransList = []
badCollectionList = []
badOrigRefList = []
noOrigRefList = []
badTransRefList = []
badMolRefList = []
badIdxRefList = []
badPclList = []
badSooList = []
badMclList = []
badSubtypeList = []
badMolMutList = []
molMutOtherNoNoteList = []

# lines that pass QC
goodLineList = []

#
# Purpose: Validate the arguments to the script.
# Returns: Nothing
# Assumes: Nothing
# Effects: sets global variable, exits if incorrect # of args
# Throws: Nothing
#
def checkArgs ():
    global inputFile

    if len(sys.argv) != 2:
        print(USAGE)
        sys.exit(1)

    inputFile = sys.argv[1]
    print('inputFile: %s' % inputFile)
    return

# end checkArgs() -------------------------------

# Purpose: create lookups, open files
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variables, exits if a file can't be opened,
#  creates files in the file system, creates connection to a database

def init ():

    # open input/output files
    openFiles()
    db.useOneConnection(1)

    #
    # create lookups
    #

    loadLookups()

    return

# end init() -------------------------------

# Purpose: load lookups for verification
# Returns: Nothing
# Assumes: 
# Effects: queries a database, modifies global variables
#
def loadLookups(): 
    global geneIdLookup, userLookup, statusLookup, typeLookup, inheritModeLookup
    global transmissionLookup, collectionLookup, referenceLookup, pclLookup
    global strainLookup, mclLookup, subtypeLookup, mutationLookup

    # Gene ID
    results = db.sql('''select accid
                from acc_accession a, mrk_marker m
                where a._mgitype_key = 2
                and a._logicaldb_key = 1
                and a.prefixPart = 'MGI:'
                and a._object_key = m._marker_key
                and m._marker_status_key in (1, 3)
                and m._organism_key = 1''', 'auto')
    for r in results:
        geneIdLookup.append(r['accid'])

    # User
    results = db.sql('''select login
                from MGI_User''', 'auto')
    for r in results:
        userLookup.append(r['login'])

    # Status
    results = db.sql('''select term
                from VOC_Term
                where _vocab_key = 37''', 'auto')
    for r in results:
        statusLookup.append(r['term'])

    # Type
    results = db.sql('''select term
                from VOC_Term
                where _vocab_key = 38''', 'auto')
    for r in results:
        typeLookup.append(r['term'])


   # Inheritance Mode
    results = db.sql('''select term
                from VOC_Term
                where _vocab_key = 35''', 'auto')
    for r in results:
        inheritModeLookup.append(r['term'])

    # Transmission
    results = db.sql('''select term
                from VOC_Term
                where _vocab_key = 61''', 'auto')
    for r in results:
        transmissionLookup.append(r['term'])

    # Collection
    results = db.sql('''select term
                from VOC_Term
                where _vocab_key = 92''', 'auto')
    for r in results:
        collectionLookup.append(r['term'])

    # Reference (JNUM)
    results = db.sql('''select accid
                from  acc_accession 
                where _mgitype_key = 1
                and _logicaldb_key = 1
                and prefixPart = 'J:'
                and preferred = 1''', 'auto')
    for r in results:
        referenceLookup.append(r['accid'])

    # Parent Cell Line name
    results = db.sql('''select cellline
                from all_cellline
                where isMutant = 0''', 'auto')
    for r in results:
        pclLookup.append(r['cellline'])

    # Strain

    results = db.sql('''select strain
                from prb_strain''', 'auto')
    for r in results:
        strainLookup.append(r['strain'])

    # Mutant Cell Line Name
    results = db.sql('''select cellline
                from all_cellline
                where isMutant = 1''', 'auto')
    for r in results:
        mclLookup.append(r['cellline'])

    # Sub Type
    results = db.sql('''select term
                from VOC_Term
                where _vocab_key = 93''', 'auto')
    for r in results:
        subtypeLookup.append(r['term'])

    # Mutation
    results = db.sql('''select term
                from VOC_Term
                where _vocab_key = 36''', 'auto')
    for r in results:
        mutationLookup.append(r['term'])
 
    return

# end loadLookups() -------------------------------

#
# Purpose: Open input and output files.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variables.
# Throws: Nothing
#
def openFiles ():
    global fpInput, fpQcRpt, fpWorkBCP, fpWarnRpt, fpDeleteRpt, fpDeleteSQL

    #
    # Open the input file
    #
    try:
        fpInput = open(inputFile, 'r')
    except:
        print('Cannot open input file: %s' % inputFile)
        sys.exit(1)

    #
    # Open QC report file
    #
    try:
        fpQcRpt = open(qcRptFile, 'w')
    except:
        print('Cannot open report file: %s' % qcRptFile)
        sys.exit(1)

    #
    # Open the warning report
    #
    try:
        fpWarnRpt = open(warnRptFile, 'w')
    except:
        print('Cannot open warning report file: %s' % warnRptFile)
        sys.exit(1)

    return

# end openFiles() -------------------------------

#
# Purpose: writes out errors to the qc report
# Returns: Nothing
# Assumes: Nothing
# Effects: writes report to the file system
# Throws: Nothing
#

def writeReport():
    global hasFatalErrors

    #
    # Now write any errors to the report
    #


    if len(missingColumnList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Lines with < 23 Columns',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(missingColumnList))

    if len(dupeLineList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Lines Duplicated',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(dupeLineList))

    if len(reqColumnList):
        hasFatalErrors = 1   
        fpQcRpt.write(CRT + CRT + str.center('Missing Data in Required Columns',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(reqColumnList))

    if len(badGeneIdList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Gene ID',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(badGeneIdList))
        
    if len(badUserList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid User Login',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(badUserList))

    if len(badTypeList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Allele Type',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(badTypeList))

    if len(badStatusList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Allele Status',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(badStatusList))

    if len(badInheritModeList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Inheritance Mode',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(badInheritModeList))

    if len(imOSNnoGenNoteList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Inheritance Mode "Other (see notes)"',60) + CRT)
        fpQcRpt.write(str.center('with no General Note',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(imOSNnoGenNoteList))

    if len(badTransList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Allele Transmission',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(badTransList))

    if len(badCollectionList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Collection',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(badCollectionList))

    if len(badOrigRefList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Original Reference',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(badOrigRefList))

    if len(noOrigRefList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Missing Original Reference',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(noOrigRefList))

    if len(badTransRefList):
        fpQcRpt.write(CRT + CRT + str.center('Invalid Transmission Reference',60)+ CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  '  + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(badTransRefList))

    if len(badMolRefList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Molecular Reference',60)+ CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(badMolRefList))

    if len(badIdxRefList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Index Reference',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(CRT.join(badIdxRefList))

    if len(badPclList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Parent Cell Line',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(CRT.join(badPclList))

    if len(badSooList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Strain of Origin',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(CRT.join(badSooList))

    if len(badMclList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Mutant Cell Line',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(CRT.join(badMclList))

    if len(badSubtypeList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Allele Subtype',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(CRT.join(badSubtypeList))

    if len(badMolMutList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Molecular Mutation',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(CRT.join(badMolMutList))

    if len(molMutOtherNoNoteList):
        hasFatalErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Molecular Mutation "Other" with no Molecular Note',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(CRT.join(molMutOtherNoNoteList))


    return

# end writeReport() -------------------------------

#
# Purpose: Close the files.
# Returns: Nothing
# Assumes: Nothing
# Effects: Modifies global variables
# Throws: Nothing
#
def closeFiles ():
    global fpInput, fpQcRpt, fpWarnRpt
    fpInput.close()
    fpQcRpt.close()
    fpWarnRpt.close()
    return

# end closeFiles) -------------------------------

def runQcChecks():
    #
    # Purpose: run all QC checks
    # Returns: Nothing
    # Assumes: Nothing
    # Effects: writes reports to the file system
    # Throws: Nothing
    #
    
    skipLine = 0
    junk = fpInput.readline() # header
    line = fpInput.readline()
    lineNum = 1
    #print('line: %s' % line)
    while line:
        lineNum += 1
        # check for dupes
        if line not in dupeLineList:
            dupeLineList.append(line)
        else:
            dupeLineList.append('%s  %s' % (lineNum, line))
            skipLine = 1
        # check that the file has at least 23 columns
        if len(str.split(line, TAB)) < 23:
            missingColumnList.append('%s  %s' % (lineNum, line))
            line = fpInput.readline()
            continue
        # get columns 1-23
        (aSym, aName, geneID, user, alleleStatus, alleleType, inheritMode, 
            transmission, collection, molNote, nomenNote, genNote, colonyNote, 
            origRef, transRef, molRef, idxRef, pcl, soo, mcl, synonym, subtype, 
            molMut) = list(map(str.strip, str.split( \
                line, TAB)))[:23]

        #
        # check required columns
        #
        if aSym == '' or aName == '' or geneID == '' or user == '' or transmission == '':
            # REPORT required fields that are empty
            reqColumnList.append('%s  %s' % (lineNum, line))
            skipLine = 1
        # 
        # verify fields that are required in the database and in file
        #
        if geneID not in geneIdLookup:
            badGeneIdList.append('%s  %s' % (lineNum, line))
            skipLine = 1

        if user not in userLookup:
            badUserList.append('%s  %s' % (lineNum, line))
            skipLine = 1            
        #
        # verify fields required in the database, but when null have defaults
        #
        if alleleStatus != '' and alleleStatus not in statusLookup:
            badStatusList.append('%s  %s' % (lineNum, line))
            skipLine = 1

        if alleleType != '': 
           if alleleType not in typeLookup:
               badTypeList.append('%s  %s' % (lineNum, line))
               skipLine = 1

        if inheritMode != '': 
            if inheritMode not in inheritModeLookup:
                badInheritModeList.append('%s  %s' % (lineNum, line))
                skipLine = 1
            elif inheritMode == 'Other (see notes)' and genNote == '':
                imOSNnoGenNoteList.append('%s  %s' % (lineNum, line))
                skipLine = 1

        if transmission != '' and transmission not in transmissionLookup:
            badTransList.append('%s  %s' % (lineNum, line))
            skipLine = 1

        if collection != '' and collection not in collectionLookup:
            badCollectionList.append('%s  %s' % (lineNum, line))
            skipLine = 1

        if origRef == '':
            noOrigRefList.append('%s  %s' % (lineNum, line))
            skipLine = 1

        elif origRef not in referenceLookup:
            badOrigRefList.append('%s  %s' % (lineNum, line))
            skipLine = 1

        if soo != '' and soo not in strainLookup:
            badSooList.append('%s  %s' % (lineNum, line))
            skipLine = 1

        #
        # Verify optional fields, some multivalued
        #

        # References (J:)
        if transRef != '' and transRef not in referenceLookup:
            badTransRefList.append('%s  %s' % (lineNum, line))

        if molRef != '' and molRef not in referenceLookup:
            badMolRefList.append('%s  %s' % (lineNum, line))
            skipLine = 1

        # can have multiple
        if idxRef != '':
            for r in str.split(idxRef, '|'):
                if r not in referenceLookup:
                    badIdxRefList.append('%s  %s' % (lineNum, line))
                    skipLine = 1

        if pcl != '' and pcl not in pclLookup:
            badPclList.append('%s  %s' % (lineNum, line))
            skipLine = 1

        # can have multiple
        if mcl != '':
            for m in str.split(mcl, '|'):
                if m not in mclLookup:
                    badMclList.append('%s  %s' % (lineNum, line))
                    skipLine = 1
 
        # can have multiple
        if subtype != '':
            for s in str.split(subtype, '|'):
                if s not in subtypeLookup:
                    badSubtypeList.append('%s  %s' % (lineNum, line))
                    skipLine = 1
            
        # can have multiple
        if molMut != '':
            # if molecular mutation = 'Other', there must be a molecular note
            for m in str.split(molMut, '|'):
                if m not in mutationLookup:
                    badMolMutList.append('%s  %s' % (lineNum, line))
                    skipLine = 1
                elif m == 'Other' and molNote == '':
                    molMutOtherNoNoteList.append('%s  %s' % (lineNum, line))
                    skipLine = 1
        if skipLine == 0:
            goodLineList.append(line)

        skipLine == 0
        line = fpInput.readline()
    return

# end runQcChecks() -------------------------------

#
# Main
#
print('checkArgs(): %s' % time.strftime("%H.%M.%S.%m.%d.%y", time.localtime(time.time())))
sys.stdout.flush()
checkArgs()

print('init(): %s' % time.strftime("%H.%M.%S.%m.%d.%y", time.localtime(time.time())))
sys.stdout.flush()
init()

print('runQcChecks(): %s' % time.strftime("%H.%M.%S.%m.%d.%y", time.localtime(time.time())))
sys.stdout.flush()
runQcChecks()

print('writeReport(): %s' % time.strftime("%H.%M.%S.%m.%d.%y", time.localtime(time.time())))
writeReport()

print('closeFiles(): %s' % time.strftime("%H.%M.%S.%m.%d.%y", time.localtime(time.time())))
sys.stdout.flush()
closeFiles()

db.useOneConnection(0)
print('done: %s' % time.strftime("%H.%M.%S.%m.%d.%y", time.localtime(time.time())))

if hasFatalErrors == 1 : 
    sys.exit(2)
else:
    sys.exit(0)

