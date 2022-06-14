
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
#      - intermediate file of QC'd alleles to create
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
#      5) Generate intermediate file of alleles to create
#      6) Close the input/output files.
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

# intermediate load ready file
loadReadyFile = os.getenv("INPUT_FILE_QC")
fpLoadReady = None

# allele types with MCLs
TAR = 'Targeted'
GT = 'Gene trapped'
EM = 'Endonuclease-mediated'
alleleTypeList = [TAR, GT, EM]

# list of line numbers written to the QC file
lineNumberSet = set([])

# MCL/PCL values
NS = 'Not Specified' # also allele type and collection default
OSN = 'Other (see notes)'

# some default values
NA = 'Not Applicable'  # inheritance mode default
RES = 'Reserved'       # allele status default

# strains with derivations in the db
one29 = '129'
one29SSvEv = '129S/SvEv'
one29P2OlaHsd = '129P2/OlaHsd'
one2955SvEvBrd = '12955/SvEvBrd'
strainList = [one29, one29SSvEv, one29P2OlaHsd, one2955SvEvBrd]

# strains in db corresponding to PCLs
# MCL NS, PCL = NS, use these PCL keys
nsOne29Key =  1098
nsOne29SSvEvKey = 40245

# MCL NS, PCL OSN, use these PCL keys
osnOne29Key = 1101
osnO29P2OlaHsdKey = 40248
osnOne2955SvEvBrdKey = 40255

# Generic 'Not Specified' PCL 
genNsPCLKey = -1

# Generic 'Other (see notes)' PCL
genOsnPCLKey = 1069

# Report file names
qcRptFile = os.environ['QC_RPT']

# 1 if any skip or warn errors in the input file
hasSkipErrors = 0
hasWarnErrors = 0

# for bcp
bcpin = '%s/bin/bcpin.csh' % os.environ['PG_DBUTILS']
server = os.environ['MGD_DBSERVER']
database = os.environ['MGD_DBNAME']

# Lookups
alleleSymbolLookup = []
geneIdLookup = {}
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
alleleInDbList = []

# allele symbols seen thus far with their line number
inputAlleleDict = {}

badGeneIdList = []
badAlleleSymbolList1 = []
badAlleleSymbolList2 = []
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

# alleles to load - they pass all QC
allelesToLoadList = []

class MutantCellLine:
    #
    # Is: data object for a mutant cell line
    # Has: if mclKey is set - use that to create association
    #      if derivationKey set - use that to create new NS mcl
    # Does: provides direct access to its attributes
    #
    def __init__(self):
        
        self.mclKeyList = []          # mutant cell line keys existing in the database
        self.derivationKey = None  # derivation key to use to create new Not Specified mcl


class Allele:
    #
    # Is: data object for a Allele
    # Has: a set of allele attributes, strings unless labeled otherwise
    # Does: provides direct access to its attributes
    #
    def __init__(self,
        aSym,           # allele symbol
        aName,          # allele name
        geneID,         # marker MGI ID 
        user,           # allele creator
        aStatus,        # allele status
        aType,          # allele type
        inheritMode,    # inheritance mode
        transmission,   # germ line transmission
        collection,     # allele collection
        molNote,        # molecular note
        nomenNote,      # allele nomenclature note
        genNote,        # general note
        colonyNote,     # pipe delim colony ID string
        origRef,        # original reference
        transRef,       # transmission reference
        molRef,         # molecular reference
        idxRefs,        # index references - multivalued '|' delimited
        synonyms,       # allele synonyms - multivalued '|' delimited
        subtypes,       # allele subtypes - multivalued '|' delimited
        molMuts,        # molecular mutations - multivalued '|' delimited
        pcl,            # parent cell line                    
        soo,            # strain of origin
        mclKeys,         # list of mcl keys with which to create allele associations 
        derivationKey): # derivation key with which to create the new Not Specified cell line

        self.aSym = aSym 
        self.aName = aName
        self.geneID = geneID
        self.user = user
        self.alleleStatus = aStatus
        self.alleleType = aType
        self.inheritMode = inheritMode
        self.transmission = transmission
        self.collection = collection 
        self.molNote = molNote
        self.nomenNote = nomenNote
        self.genNote = genNote
        self.colonyNote = colonyNote
        self.origRef = origRef
        self.transRef = transRef
        self.molRef = molRef
        self.idxRefs = idxRefs 
        self.synonyms = synonyms
        self.subtypes = subtypes
        self.molMuts = molMuts
        self.pcl = pcl
        self.soo = soo
        self.mclKeys = mclKeys
        self.derivationKey = derivationKey

    def toString(this):
        return '%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s' % (this.aSym, this.aName, this.geneID, this.user, this.alleleStatus, this.alleleType, this.inheritMode, this.transmission, this.collection, this.molNote, this.nomenNote, this.genNote, this.colonyNote, this.origRef, this.transRef, this.molRef, this.idxRefs, this.synonyms, this.subtypes, this.molMuts, this.pcl, this.soo, this.mclKeys, this.derivationKey)

    def toLoad(this):
        return '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (this.aSym, TAB, this.aName, TAB, this.geneID, TAB, this.user, TAB, this.alleleStatus, TAB, this.alleleType, TAB, this.inheritMode, TAB, this.transmission, TAB, this.collection, TAB, this.molNote, TAB, this.nomenNote, TAB, this.genNote, TAB, this.colonyNote, TAB, this.origRef, TAB, this.transRef, TAB, this.molRef, TAB, this.idxRefs, TAB, this.synonyms, TAB, this.subtypes, TAB, this.molMuts, TAB, this.pcl, TAB, this.soo, TAB, this.mclKeys, TAB, this.derivationKey, CRT)

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
    results = db.sql('''select a.accid, m.symbol
                from acc_accession a, mrk_marker m
                where a._mgitype_key = 2
                and a._logicaldb_key = 1
                and a.prefixPart = 'MGI:'
                and a._object_key = m._marker_key
                and m._marker_status_key in (1, 3)
                and m._organism_key = 1''', 'auto')
    for r in results:
        geneIdLookup[r['accid']] = r['symbol']

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
    global fpInput, fpLoadReady, fpQcRpt

    #
    # Open the input file
    #
# encoding='utf-8' no
# encoding=u'utf-8' no

    try:
        fpInput = open(inputFile, 'r', encoding='utf-8', errors='replace')
    except:
        print('Cannot open input file: %s' % inputFile)
        sys.exit(1)

    #
    # Open load ready input file
    #
    try:
        fpLoadReady = open(loadReadyFile, 'w')
    except:
        print('Cannot open load ready file: %s' % loadReadyFile)
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
    global hasSkipErrors, hasWarnErrors, lineNumberSet

    #
    # Now write any errors to the report
    #
    fpQcRpt.write( str.center('Warning QC - these will be loaded',80) + CRT)

    if len(alleleInDbList):
        hasWarnErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Allele Symbols already in the DB (case sensitive)',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(alleleInDbList))
        fpQcRpt.write(CRT + 'Total: %s' % len(alleleInDbList))

    # iterate thru the dictionary of all symbols/line numbers in input
    found = 0
    dupeSymCount = 0
    for a in inputAlleleDict: # contains all alleles in input, with a list of their line numbers
        if len(inputAlleleDict[a]) > 1:  # we only want to report if the symbol found more than once
            if found == 0: # print the report section header if we found
                fpQcRpt.write(CRT + CRT + str.center('Allele Symbols duplicated in the input file (case sensitive)',60) + CRT)
                fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
                fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
                found = 1
            hasWarnErrors = 1
            fpQcRpt.write('%s    ' % (a))
            fpQcRpt.write(', '.join(inputAlleleDict[a]))
            fpQcRpt.write(CRT)
            dupeSymCount +=1
    if found > 0:
        fpQcRpt.write(CRT + 'Total: %s' % dupeSymCount)

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

    if len(badAlleleSymbolList1):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Marker symbol not in Allele symbol and Allele not Transgenic',60) + CRT)
        fpQcRpt.write('%-12s  %-20s  %-20s%s' % ('Line#','MSymbol', 'Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + 20*'-' + CRT)
        fpQcRpt.write(''.join(badAlleleSymbolList1))
        fpQcRpt.write(CRT + 'Total: %s' % len(badAlleleSymbolList1))

    if len(badAlleleSymbolList2):
        hasSkipErrors = 1
        fpQcRpt.write(CRT + CRT + str.center('Allele symbol must either have both < and > or neither',60) + CRT)
        fpQcRpt.write('%-12s  %-20s%s' % ('Line#','Line', CRT))
        fpQcRpt.write(12*'-' + '  ' + 20*'-' + CRT)
        fpQcRpt.write(''.join(badAlleleSymbolList2))
        fpQcRpt.write(CRT + 'Total: %s' % len(badAlleleSymbolList2))

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

    # for development - helps to determine which lines in test file are being reported
    #fpQcRpt.write(CRT + CRT + 'sorted list of line numbers reported: ' + CRT)
    #sortedList =  list(lineNumberSet)
    #sortedList.sort()
    #s = [str(i) for i in sortedList]
    #fpQcRpt.write(', '.join(s))
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
    global fpInput, fpLoadReady, fpQcRpt
    fpInput.close()
    fpLoadReady.close()
    fpQcRpt.close()

    return

# end closeFiles) -------------------------------

    #
    # Purpose: run all QC checks
    # Returns: Nothing
    # Assumes: file descriptors have been initialized
    # Effects: writes reports and the load ready file to file system
    # Throws: Nothing
    #

def runQcChecks():
    
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
        #print('lineNum: %s %s' % (lineNum, line))
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
            origRef, transRef, molRef, idxRefs, pcl, soo, mcls, synonyms, 
            subtypes, molMuts) = list(map(str.strip, str.split( \
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
        if alleleType in [TAR, GT] and (not mcls or not pcl):
            tarGtMissingMclPclList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)

        # if allele type is EM then both MCL and PCL must be specified
        # matrix rows 4/5b
        if alleleType == EM:
            if (mcls and not pcl) or (not mcls and pcl): 
                emMissingMclPclList.append('%s  %s' % (lineNum, line))
                skipLine = 1
                lineNumberSet.add(lineNum)

        # Non TAR/GT/EM alleles should have neither MCL or PCL specified
        # matrix rows 7/8
        if not badAlleleType and alleleType not in alleleTypeList and (mcls or pcl):
            nonTARGTEMwithMclPclList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)
        # 
        # verify fields that are required in the database and in file
        #

        # if neither < or > OK, if both < and > OK, otherwise report
        if (aSym.find('<') == -1 and aSym.find('>') == -1) or (aSym.find('<') != -1 and aSym.find('>') != -1):
            pass
        else:
            badAlleleSymbolList2.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)

        if aSym in alleleSymbolLookup:
            alleleInDbList.append('%s  %s' % (lineNum, line))
            lineNumberSet.add(lineNum)

        if aSym not in inputAlleleDict:
            inputAlleleDict[aSym] = []
        inputAlleleDict[aSym].append(str(lineNum))

        if geneID not in geneIdLookup:
            badGeneIdList.append('%s  %s' % (lineNum, line))
            skipLine = 1
            lineNumberSet.add(lineNum)

        # if the marker symbol is not part of the allele symbol and allele type is not transgenic
        if aSym.find(geneIdLookup[geneID]) == -1 and alleleType != 'Transgenic':
            # report and skip
            badAlleleSymbolList1.append('%s  %s  %s\n' % (lineNum, geneIdLookup[geneID], line))
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
        if idxRefs != '':
            for r in str.split(idxRefs, '|'):
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
        if mcls != '':
            for m in str.split(mcls, '|'):
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
        # can have multiple
        if subtypes != '':
            for s in str.split(subtypes, '|'):
                if s not in subtypeLookup:
                    badSubtypeList.append('%s  %s' % (lineNum, line))
                    skipLine = 1
                    lineNumberSet.add(lineNum) 
        # can have multiple
        if molMuts != '':
            # if molecular mutation = 'Other', there must be a molecular note
            for m in str.split(molMuts, '|'):
                if m not in mutationLookup:
                    badMolMutList.append('%s  %s' % (lineNum, line))
                    skipLine = 1
                    lineNumberSet.add(lineNum)
                elif m == 'Other' and molNote == '':
                    molMutOtherNoNoteList.append('%s  %s' % (lineNum, line))
                    skipLine = 1
                    lineNumberSet.add(lineNum)
        if skipLine == 0:
            # A list of MutantCellLine objects
            resolvedMcls = []

            # these will remain blank if not TAR/GT or EM allele type
            mclKeys = ''
            derivationKey = ''

            # We have a TAR/GT/EM allele that passes QC thus far AND mcls and
            # pcl are specified (EM can have null mcl/pcl in which case there 
            # will be no mcl association/creation)
            # now find the mcl or correct derivation, to create NS MCL, 
            # using the MCL/PCL/SOO rules
            if alleleType in alleleTypeList and mcls and pcl:
                # attempt to resolve the mcls for this allele
                # if resolvedMcls empty, we could not resolve so skip this allele
                resolvedMcls = qcMCL(aSym, mcls, pcl, soo, alleleType, line, lineNum) 
                if not resolvedMcls:
                    skipLine = 1
                else:
                    for mObject in resolvedMcls:
                        if mObject.mclKeyList:
                            mclKeys = '|'.join(mObject.mclKeyList)
                        elif mObject.derivationKey:
                            derivationKey = mObject.derivationKey
        if skipLine == 0:
            goodLineList.append(line)
            if alleleStatus == '':
                alleleStatus = RES
            if alleleType == '':
                alleleType = NS
            if inheritMode == '':
                inheritMode = NA
            if collection == '':
                collection = NS
            
            
            alleleToLoad = Allele(aSym, aName, geneID, user, alleleStatus, alleleType, inheritMode, transmission, collection, molNote, nomenNote, genNote, colonyNote, origRef, transRef, molRef, idxRefs, synonyms, subtypes, molMuts, pcl, soo, mclKeys, derivationKey)
            allelesToLoadList.append(alleleToLoad)
            #print('%s %s' % (lineNum, alleleToLoad.toString()))
        skipLine = 0

        line = fpInput.readline()
    return

# end runQcChecks() -------------------------------

#
# Purpose: QC the MCL and a) find MCL in database to associated with the 
#       allele or find the derivation in the database with which to create
#       a new Not Specified MCL
# Returns:  a list of MutantCellLine objects
# Assumes: Nothing
# Effects:  queries a database
# Throws: Nothing
#

def qcMCL(aSym, mcls, pcl, soo, alleleType, line, lineNum):

    # a list of MutantCellLine objects (see class)
    resolvedMclList = []

    #print(CRT + CRT + 'In qcMCL')
    #print('in qcMCL  lineNum: %s allele symbol: %s, mcls: %s pcl: %s soo: %s alleleType: %s ' % (lineNum, aSym, mcls, pcl, soo, alleleType))

    for m in str.split(mcls, '|'):
        #print('m: %s' % m)
        if m != NS: # rows 10-12 in the matrix
            #print('mcl != NS: lineNum: %s allele symbol: %s, mcls: %s pcl: %s soo: %s alleleType: %s' % (lineNum, aSym, mcls, pcl, soo, alleleType))
            # lookup PCL for MCL in ALL_CellLine_Derivation_view
            # if same as incoming PCL and incoming strain, QC passes
            sql = '''select c._cellline_key, v.parentcellline, v.parentcelllinestrain
                from all_cellline c, all_cellLine_derivation_view v, voc_term t
                where c.isMutant = 1
                and c.cellline = '%s'
                and v._derivationtype_key = t._term_key
                and c._derivation_key = v._derivation_key''' % (m)

            #print(sql)
            results = db.sql(sql, 'auto')
            
            if len(results) != 1:
                print ('result != 1 for mcl: %s  %s' % (m, results))
            else:
                dbPcl = results[0]['parentcellline']
                dbStrain =  results[0]['parentcelllinestrain']
                mclKey = results[0]['_cellline_key']
                # if the incoming pcl does not match the mcl pcl in the database
                # report and skip
                if pcl != dbPcl:
                    pclNEdbPclList.append('%s  %s' % (lineNum, line))
                    lineNumberSet.add(lineNum)
                    #print('pcl != dbPcl')

                # if the incoming soo does not match the pcl strain in the database
                # report and skip
                elif soo != dbStrain:
                    sooNEdbStrainList.append('%s  %s' % (lineNum, line))
                    lineNumberSet.add(lineNum)
                    #print('soo != dbStrain')
                
                # otherwise use the incoming named mcl
                else:
                    #print('mcl != NS - incoming mcl used')
                    mclResolved = MutantCellLine()
                    mclResolved.mclKeyList.append(str(mclKey))
                    resolvedMclList.append(mclResolved)
                
        else: # m == NS
            if pcl not in (NS, OSN):
                #print('mcl=NS, pcl not in (NS, OSN): lineNum: %s allele symbol: %s, mcls: %s pcl: %s soo: %s alleleType: %s' % (lineNum, aSym, mcls, pcl, soo, alleleType))
                # find the PCL
                results = db.sql('''select c._cellline_key as _parentcellline_key, c.cellline as parentcellline, c.celllinestrain as parentcelllinestrain
                    from all_cellline_view c
                    where c.isMutant = 0
                    and cellline = '%s' ''' % pcl, 'auto')
                #print('length of results: %s' % len(results))
                # check that the pcl strain in db same as incoming soo
                if soo != results[0]['parentcelllinestrain']:
                    sooNEdbStrainList.append('%s  %s' % (lineNum, line))
                    lineNumberSet.add(lineNum)
                    #print('soo != dbStrain')
                else:
                    pclKey = results[0]['_parentcellline_key']
                    # find the derivation to use
                    sql = '''select v.name, v._derivation_key, 
                            v._parentcellline_key, v.parentcellline, v.parentcelllinestrain 
                        from all_cellline_derivation_view v, voc_term t
                        where v._parentcellline_key = %s
                         and v.creator = '%s'
                        and v._derivationtype_key = t._term_key
                        and t.term = '%s' ''' % (pclKey, NS, alleleType)
                    #print(sql + '\n')
                    results = db.sql(sql, 'auto')
                    #print(results)
                    #print(CRT + CRT)   
                    mclToCreate = MutantCellLine()
                    mclToCreate.derivationKey = results[0]['_derivation_key']
                    resolvedMclList.append(mclToCreate)
            elif pcl == NS:
                # no checking needed here - just have to determine the correct Derivation to create the 
                # new NS MCL with
                #print ('mcl NS, pcl NS: lineNum: %s allele symbol: %s, mcls: %s pcl: %s soo: %s alleleType: %s ' % (lineNum, aSym, mcls, pcl, soo, alleleType))
                # find the pcl
                pclKeyToUse = 0

                if soo == one29:
                    pclKeyToUse = osnOne29Key
                elif soo == one29SSvEv:
                    pclKeyToUse = nsOne29SSvEvKey
                else:
                    pclKeyToUse = genNsPCLKey

                sql = '''select v.name, v._derivation_key,
                            v._parentcellline_key, v.parentcellline, v.parentcelllinestrain
                        from all_cellline_derivation_view v, voc_term t
                        where v._parentcellline_key = '%s'
                        and v.creator = '%s'
                        and v._derivationtype_key = t._term_key
                        and t.term = '%s' ''' % (pclKeyToUse, NS, alleleType)
                #print(sql + '\n')
                results = db.sql(sql, 'auto')
                #print(results)
                #print(CRT + CRT)
                mclToCreate = MutantCellLine()
                mclToCreate.derivationKey = results[0]['_derivation_key']
                resolvedMclList.append(mclToCreate)

# strains with derivations in the db
            elif pcl == OSN:
                # no checking needed here - just have to determine the correct Derivation to create the
                # new NS MCL with  
                #print ('mcl NS, pcl OSN: lineNum: %s allele symbol: %s, mcls: %s pcl: %s soo: %s alleleType: %s ' % (lineNum, aSym, mcls, pcl, soo, alleleType))
                # find the pcl
                pclKeyToUse = 0

                if soo == one29:
                    pclKeyToUse = osnOne29Key
                elif soo == one29P2OlaHsd:
                    pclKeyToUse = osnO29P2OlaHsdKey
                elif soo == one2955SvEvBrd:
                    pclKeyToUse = osnOne2955SvEvBrdKey
                else:
                    pclKeyToUse = genOsnPCLKey

                # now find the derivation
                sql = '''select v.name, v._derivation_key,
                            v._parentcellline_key, v.parentcellline, v.parentcelllinestrain
                        from all_cellline_derivation_view v, voc_term t
                        where v._parentcellline_key = '%s'
                        and v.creator = '%s'
                        and v._derivationtype_key = t._term_key
                        and t.term = '%s' ''' % (pclKeyToUse, NS, alleleType)
                #print(sql + '\n')
                results = db.sql(sql, 'auto')
                #print(results)
                #print(CRT + CRT)
                mclToCreate = MutantCellLine()
                mclToCreate.derivationKey = results[0]['_derivation_key']
                # wait - are we creating a mcl for this?
                resolvedMclList.append(mclToCreate)
    return resolvedMclList

# end qcMCL() -------------------------------

def writeLoadReadyFile():
    for a in allelesToLoadList:
        fpLoadReady.write(a.toLoad())

    return

# end writeLoadReadyFile() -------------------------------

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

print('writeLoadReadyFile(): %s' % time.strftime("%H.%M.%S.%m.%d.%y", time.localtime(time.time())))
writeLoadReadyFile()

print('closeFiles(): %s' % time.strftime("%H.%M.%S.%m.%d.%y", time.localtime(time.time())))
sys.stdout.flush()
closeFiles()

db.useOneConnection(0)
print('done: %s' % time.strftime("%H.%M.%S.%m.%d.%y", time.localtime(time.time())))
if hasSkipErrors and hasWarnErrors:
    sys.exit(2)
elif hasSkipErrors: 
    sys.exit(3)
elif hasWarnErrors:
    sys.exit(4)
else:
    sys.exit(0)

