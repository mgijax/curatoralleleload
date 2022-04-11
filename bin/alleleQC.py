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
#      - temp table BCP file (${MGI_ID_BCP})
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#      2:  Report/Skip and/or Warning QC errors detected and written to report
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
import Set

#
#  CONSTANTS
#
TAB = '\t'
CRT = '\n'

USAGE = 'Usage: alleleQC.py  inputFile'

#
#  GLOBALS
#

# allele types with MCLs
TAR = 'Targeted'
GT = 'Gene trapped'
EM = 'Endonuclease-mediated'
alleleTypeList = [TAR, GT, EM]

# list of line numbers written to the QC file
lineNumberSet = set([])

# MCL/PCL values
NS = 'Not Specified'
OSN = 'Other (see notes)'

# strains with derivations in the db
one29 = '129'
one29SSvEv = '129S/SvEv'
one29P2OlaHsd = '129P2/OlaHsd'
one2955SvEvBrd = '12955/SvEvBrd'
strainList = [one29, one29SSvEv, one29P2OlaHsd, one2955SvEvBrd]

# strains in db corresponding to PCLs
# MCL NS, PCL = NS, use these PCL keys
nsOne29 =  1098
nsOne29SSvEv = 40245

# MCL NS, PCL OSN, use these PCL keys
osnOne29 = 1101
osnO29P2OlaHsd = 40248
osnOne2955SvEvBrd = 40255

# Report file names
qcRptFile = os.environ['QC_RPT']

# 1 if any QC errors in the input file
hasSkipErrors = 0
hasWarnErrors = 0

# for bcp
bcpin = '%s/bin/bcpin.csh' % os.environ['PG_DBUTILS']
server = os.environ['MGD_DBSERVER']
database = os.environ['MGD_DBNAME']

# Lookups
alleleSymbolLookup = []
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
tarGtMissingMclPclList = []
emMissingMclPclList = []
nonTARGTEMwithMclPclList = []
dupeAlleleList = []
badGeneIdList = []
badUserList = []
badStatusList = []
badTypeList = []
badInheritModeList = []
# inheritance mode is 'Other (see notes)' and there is not general note
imOSNnoGenNoteList = [] 
badTransList = []
badCollectionList = []
noOrigRefList = []
badOrigRefList = []
badTransRefList = []
badMolRefList = []
badIdxRefList = []
badPclList = []
osnPclNoGenNoteList = []
badSooList = []
badMclList = []
mismatchedGeneIDList = []
badSubtypeList = []
badMolMutList = []
molMutOtherNoNoteList = []

pclNEdbPclList = []
sooNEdbStrainList = []

# lines seen in the input file
distinctLineList = []

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
    global alleleSymbolLookup, geneIdLookup, userLookup, statusLookup, typeLookup, inheritModeLookup
    global transmissionLookup, collectionLookup, referenceLookup, pclLookup
    global strainLookup, mclLookup, subtypeLookup, mutationLookup

    # Allele Symbol
    results = db.sql('''select distinct symbol
                from all_allele''', 'auto')

    for r in results:
        alleleSymbolLookup.append(r['symbol'])

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
    global fpInput, fpQcRpt, fpWorkBCP, fpDeleteRpt, fpDeleteSQL

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
    global hasSkipErrors, lineNumberSet

    #
    # Now write any errors to the report
    #
    fpQcRpt.write( str.center('Warning QC - these will be loaded',80) + CRT)

    if len(dupeAlleleList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Allele Symbols already in the DB (case sensitive)',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(dupeAlleleList))
        fpQcRpt.write(CRT + 'Total: %s' % len(dupeAlleleList))
    
    fpQcRpt.write(CRT + CRT)
    fpQcRpt.write( str.center('Report/Skip QC - these will be reported and skipped',80) + CRT)

    if len(dupeLineList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Lines Duplicated',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(dupeLineList))
        fpQcRpt.write(CRT + 'Total: %s' % len(dupeLineList))

    if len(missingColumnList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Lines with < 23 Columns',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(CRT.join(missingColumnList))
        fpQcRpt.write(CRT + 'Total: %s' % len(missingColumnList))

    if len(reqColumnList):
        hasSkipErrors = 1   
        fpQcRpt.write(CRT + CRT + str.center('Missing Data in Required Columns',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(reqColumnList))
        fpQcRpt.write(CRT + 'Total: %s' % len(reqColumnList))

    if len(tarGtMissingMclPclList):
        fpQcRpt.write(CRT + CRT + str.center('TAR/GT Allele with missing MCL or PCL', 60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(tarGtMissingMclPclList))
        fpQcRpt.write(CRT + 'Total: %s' % len(tarGtMissingMclPclList))

    if len(emMissingMclPclList):
        fpQcRpt.write(CRT + CRT + str.center('EM Allele with missing MCL or PCL', 60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(emMissingMclPclList))
        fpQcRpt.write(CRT + 'Total: %s' % len(emMissingMclPclList))

    if len(nonTARGTEMwithMclPclList):
        fpQcRpt.write(CRT + CRT + str.center('Non TAR/GT/EM Allele with specified MCL and/or PCL', 60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(nonTARGTEMwithMclPclList))
        fpQcRpt.write(CRT + 'Total: %s' % len(nonTARGTEMwithMclPclList))

    if len(badGeneIdList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Gene ID',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(badGeneIdList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badGeneIdList))
        
    if len(badUserList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid User Login',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(badUserList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badUserList))

    if len(badStatusList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Allele Status',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(badStatusList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badStatusList))

    if len(badTypeList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Allele Type',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(badTypeList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badTypeList))

    if len(badInheritModeList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Inheritance Mode',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(badInheritModeList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badInheritModeList))

    if len(imOSNnoGenNoteList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Inheritance Mode "Other (see notes)"',60) + CRT)
        fpQcRpt.write(str.center('with no General Note',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(imOSNnoGenNoteList))
        fpQcRpt.write(CRT + 'Total: %s' % len(imOSNnoGenNoteList))

    if len(badTransList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Allele Transmission',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(badTransList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badTransList))

    if len(badCollectionList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Collection',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(badCollectionList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badCollectionList))

    if len(noOrigRefList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Missing Original Reference',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(noOrigRefList))
        fpQcRpt.write(CRT + 'Total: %s' % len(noOrigRefList))

    if len(badOrigRefList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Original Reference',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(badOrigRefList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badOrigRefList))

    if len(badTransRefList):
        fpQcRpt.write(CRT + CRT + str.center('Invalid Transmission Reference',60)+ CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  '  + 20*'-' + CRT)
        fpQcRpt.write(''.join(badTransRefList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badTransRefList))

    if len(badMolRefList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Molecular Reference',60)+ CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(badMolRefList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badMolRefList))

    if len(badIdxRefList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Index Reference',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(''.join(badIdxRefList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badIdxRefList))

    if len(badPclList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Parent Cell Line',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(''.join(badPclList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badPclList))

    if len(osnPclNoGenNoteList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('PCL Other (see notes) with no General Note',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(''.join(osnPclNoGenNoteList))
        fpQcRpt.write(CRT + 'Total: %s' % len(osnPclNoGenNoteList))

    if len(badSooList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Strain of Origin',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(''.join(badSooList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badSooList))

    if len(badMclList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Mutant Cell Line',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(''.join(badMclList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badMclList))

    if len(mismatchedGeneIDList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center("MCL marker doesn't match input marker",60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(''.join(mismatchedGeneIDList))
        fpQcRpt.write(CRT + 'Total: %s' % len(mismatchedGeneIDList))

    if len(badSubtypeList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Allele Subtype',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(''.join(badSubtypeList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badSubtypeList))

    if len(badMolMutList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Invalid Molecular Mutation',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(''.join(badMolMutList))
        fpQcRpt.write(CRT + 'Total: %s' % len(badMolMutList))

    if len(molMutOtherNoNoteList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Molecular Mutation "Other" with no Molecular Note',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(''.join(molMutOtherNoNoteList))
        fpQcRpt.write(CRT + 'Total: %s' % len(molMutOtherNoNoteList))

    if len(pclNEdbPclList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Specified MCL where input PCL != DB PCL',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(''.join(pclNEdbPclList))
        fpQcRpt.write(CRT + 'Total: %s' % len(pclNEdbPclList))

    if len(sooNEdbStrainList):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Specified MCL where input SOO != DB PCL Strain',60)+ CRT)
        fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
        fpQcRpt.write(''.join(sooNEdbStrainList))
        fpQcRpt.write(CRT + 'Total: %s' % len(sooNEdbStrainList))

    fpQcRpt.write(CRT + CRT + 'sorted list of line numbers reported: ' + CRT)
    sortedList =  list(lineNumberSet)
    sortedList.sort()
    s = [str(i) for i in sortedList]
    fpQcRpt.write(', '.join(s))
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
    global fpInput, fpQcRpt
    fpInput.close()
    fpQcRpt.close()
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
    global lineNumberSet
 
    skipLine = 0
    junk = fpInput.readline() # header
    line = fpInput.readline()
    lineNum = 1
    while line:

        # flag so we don't report lines with bad allele type n the Non TAR/GT/EM Allele
        # with specified MCL and/or PCL section
        badAlleleType = 0
        lineNum += 1
        #print('line: %s %s' % (lineNum, line))
        # check for dupes
        if line not in distinctLineList:
            distinctLineList.append(line)
        else:
            dupeLineList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)
        # check that the file has at least 23 columns
        if len(str.split(line, TAB)) < 23:
            missingColumnList.append('%s  %s' % (lineNum, line))
            line = fpInput.readline()
            lineNumberSet.add(lineNum)
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
        if aSym == '' or aName == '' or geneID == '' or user == '' or transmission == '' or soo == '':
            # REPORT required fields that are empty
            reqColumnList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)

        # alleleType has default when empty
        if alleleType != '' and alleleType not in typeLookup:
              badTypeList.append('%s  %s' % (lineNum, line))
              skipLine = 1
              lineNumberSet.add(lineNum)
              badAlleleType = 1 

        # if allele type is TAR/GT and MCL is null or PCL is null
        # matrix rows 4/5a
        if alleleType in [TAR, GT] and (not mcl or not pcl):
            tarGtMissingMclPclList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)

        # if allele type is EM then both MCL and PCL must be specified
        # matrix rows 4/5b
        if alleleType == EM:
            if (mcl and not pcl) or (not mcl and pcl): 
                emMissingMclPclList.append('%s  %s' % (lineNum, line))
                skipLine = 1
                lineNumberSet.add(lineNum)

        # Non TAR/GT/EM alleles should have neither MCL or PCL specified
        # matrix rows 7/8
        if not badAlleleType and alleleType not in alleleTypeList and (mcl or pcl):
            nonTARGTEMwithMclPclList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)
        # 
        # verify fields that are required in the database and in file
        #
        if aSym in alleleSymbolLookup:
            dupeAlleleList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)

        if geneID not in geneIdLookup:
            badGeneIdList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)
        if user not in userLookup:
            badUserList.append('%s  %s' % (lineNum, line))
            skipLine = 1            
            lineNumberSet.add(lineNum)
        #
        # verify fields required in the database, but when null have defaults
        #
        if alleleStatus != '' and  alleleStatus not in statusLookup:
            badStatusList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)
        if inheritMode != '': 
            if inheritMode not in inheritModeLookup:
                badInheritModeList.append('%s  %s' % (lineNum, line))
                skipLine = 1
                lineNumberSet.add(lineNum)
            elif inheritMode == OSN and genNote == '':
                imOSNnoGenNoteList.append('%s  %s' % (lineNum, line))
                skipLine = 1
                lineNumberSet.add(lineNum)
        if transmission != '' and transmission not in transmissionLookup:
            badTransList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)

        if collection != '' and collection not in collectionLookup:
            badCollectionList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)
        if origRef == '':
            noOrigRefList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)
        elif origRef not in referenceLookup:
            badOrigRefList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)
        if soo != '' and soo not in strainLookup:
            badSooList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)
        #
        # Verify optional fields, some multivalued
        #

        # References (J:)
        if transRef != '' and transRef not in referenceLookup:
            badTransRefList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)
        if molRef != '' and molRef not in referenceLookup:
            badMolRefList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)
        # can have multiple
        if idxRef != '':
            for r in str.split(idxRef, '|'):
                if r not in referenceLookup:
                    badIdxRefList.append('%s  %s' % (lineNum, line))
                    skipLine = 1
                    lineNumberSet.add(lineNum)
        if pcl != '' and pcl not in pclLookup:
            badPclList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)
        if pcl!= '' and pcl == OSN and genNote == '':
            osnPclNoGenNoteList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)
        # can have multiple
        if mcl != '':
            for m in str.split(mcl, '|'):
                if m not in mclLookup:
                    #print('bad mcl: %s' % m)
                    badMclList.append('%s  %s' % (lineNum, line))
                    skipLine = 1
                    lineNumberSet.add(lineNum) 
                # if MCL in the database, lookup it's marker ID in the db, report
                # if different than the incoming marker ID
                elif m != NS:
                    sql = '''select a.accid
                        from all_allele_cellLine_view v, acc_accession a, all_allele aa
                        where v.isMutant = 1
                        and v.cellline = '%s'
                        and v._allele_key = aa._allele_key
                        and aa._marker_key = a._object_key
                        and a._mgitype_key = 2
                        and a.preferred = 1
                        and a._logicaldb_key = 1 ''' % m
                    #print(sql)
                    results = db.sql(sql, 'auto')
                    #print(results)
                    if len(results) < 1 or len(results) > 1:
                        print('MCL is not NS and marker id lookup has no results or too many results')
                        #print(results)
                    else:
                        dbGeneID = results[0]['accid']
                        if geneID != dbGeneID:
                            mismatchedGeneIDList.append('%s  %s' % (lineNum, line))
                            skipLine = 1
                            lineNumberSet.add(lineNum)
                            #print('yes marker mismatch')
                        #else:
                        #    print('no marker mismatch')
        # can have multiple
        if subtype != '':
            for s in str.split(subtype, '|'):
                if s not in subtypeLookup:
                    badSubtypeList.append('%s  %s' % (lineNum, line))
                    skipLine = 1
                    lineNumberSet.add(lineNum) 
        # can have multiple
        if molMut != '':
            # if molecular mutation = 'Other', there must be a molecular note
            for m in str.split(molMut, '|'):
                if m not in mutationLookup:
                    badMolMutList.append('%s  %s' % (lineNum, line))
                    skipLine = 1
                    lineNumberSet.add(lineNum)
                elif m == 'Other' and molNote == '':
                    molMutOtherNoNoteList.append('%s  %s' % (lineNum, line))
                    skipLine = 1
                    lineNumberSet.add(lineNum)
        if skipLine == 0:
            # We have a TAR/GT/EM allele that passes QC thus far - now check
            # MCL/PCL/SOO rules
            if alleleType in alleleTypeList:
                skipLine = qcMCL(aSym, mcl, pcl, soo, alleleType, genNote, line, lineNum)       

        if skipLine == 0:
            goodLineList.append(line)
        skipLine = 0

        line = fpInput.readline()
    return

# end runQcChecks() -------------------------------
# strains with derivations in the db
#one29 = '129'
#one29SSvEv = '129S/SvEv'
#one29P2OlaHsd = '129P2/OlaHsd'
#one2955SvEvBrd = '12955/SvEvBrd'
#strainList = [one29, one29SSvEv, one29P2OlaHsd, one2955SvEvBrd]

# strains in db corresponding to PCLs
# MCL NS, PCL = NS, use these PCL keys
#nsOne29 =  1098
#nsOne29SSvEv = 40245

# MCL NS, PCL OSN, use these PCL keys
#osnOne29 = 1101
#osnO29P2OlaHsd = 40248
#osnOne2955SvEvBrd = 40255

def qcMCL(aSym, mcl, pcl, soo, alleleType, genNote, line, lineNum):
    #global 
    skipLine = 0
    print(CRT + CRT + 'In qcMCL')
    print('lineNum: %s allele symbol: %s, mcl: %s pcl: %s soo: %s alleleType: %s ' % (lineNum, aSym, mcl, pcl, soo, alleleType))
    for m in str.split(mcl, '|'):
        print('m: %s' % m)
        if m != NS: # rows 10-12 in the matrix
            print('row 10 checks: m != NS')
            # lookup PCL for MCL in ALL_CellLine_Derivation_view
            # if same as incoming PCL and incoming strain, QC passes
            sql = '''select v.parentcellline, v.parentcelllinestrain
                from all_cellline c, all_cellLine_derivation_view v, voc_term t
                where c.isMutant = 1
                and c.cellline = '%s'
                and v._derivationtype_key = t._term_key
                and c._derivation_key = v._derivation_key''' % (m)

            print(sql)
            results = db.sql(sql, 'auto')
            
            if len(results) != 1:
                print ('result != 1 for mcl: %s  %s' % (m, results))
            else:
                dbPcl = results[0]['parentcellline']
                dbStrain =  results[0]['parentcelllinestrain']
                if pcl != dbPcl:
                    pclNEdbPclList.append('%s  %s' % (lineNum, line))
                    lineNumberSet.add(lineNum)
                    skipLine = 1
                    print('pcl != dbPcl')
                elif soo != dbStrain:
                    sooNEdbStrainList.append('%s  %s' % (lineNum, line))
                    lineNumberSet.add(lineNum)
                    skipLine = 1
                    print('soo != dbStrain')
                print(results)
                print(CRT + CRT)
                #else:
                    # determine the correct PCL to create the new NS MCL with
                #    if soo == 

        else: # m == NS
            if pcl not in (NS, OSN):
                print('mcl=NS, pcl not in (NS, OSN)')
                # find the PCL
                results = db.sql('''select c.cellline as parentcellline, c.celllinestrain parentcelllinestrain
                    from all_cellline_view c
                    where c.isMutant = 0
                    and cellline = '%s' ''' % pcl, 'auto')
                print('length of results: %s' % len(results))
                # check that the pcl strain in db same as incoming soo
                if soo != results[0]['parentcelllinestrain']:
                    sooNEdbStrainList.append('%s  %s' % (lineNum, line))
                    lineNumberSet.add(lineNum)
                    skipLine = 1
                    print('soo != dbStrain')
                #else:
                    # find the MCL to use
                print(results)
                print(CRT + CRT)   
            #else: # pcl in (NS, OSN)
                # no checking needed here - just have to determine the correct PCL to create the 
                # new NS MCL with

    return skipLine

# end qcMCL() -------------------------------

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

if hasSkipErrors == 1 : 
    sys.exit(2)
else:
    sys.exit(0)

