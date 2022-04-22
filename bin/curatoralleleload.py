#
# curatoralleleload.py
###############################################################################
#
# Purpose:
#
#	An curator driven allele load to load alleles and their accids, notes,
#           references, synonyms, mutations, subtypes, mcl's and  cellline associations
#
# Usage:
#	curatoralleleload.py 
#
# Envvars:
#
# Inputs:
#
#	A tab-delimited file in the format:
#
#       field 1: Allele Symbol
#       field 2: Allele Name
#	field 3: MGI Gene ID
#       field 4: Created By
#       field 5: Allele Status
#       field 6: Allele Type
#       field 7: Inheritance Mode
#       field 8: Transmission
#       field 9: Allele Collection
#       field 10: Molecular Note 
#       field 11: Nomenclature Note
#       field 12: General Note
#       field 13: Colony ID Note
#       field 14: Original Reference
#       field 15: Transmission Reference
#       field 16: Molecular Reference
#       field 17: Index References - multivalued '|' delimited
#       field 18: Synonyms - General. multivalued '|' delimited
#       field 19: Allele Subtypes - multivalued '|' delimited
#       field 20: Molecular Mutations, multivalued '|' delimited
#       field 21: Parent Cell Line
#       field 22: Strain of Origin Name
#       field 23: Mutant Cell Lines - multivalued '|' delimited
#
# Outputs:
#
#       BCP files:
#
#       ALL_Allele.bcp                  master Allele records
#       ACC_Accession.bcp               Accession records
#       MGI_Note                        allele notes (all types)
#       MGI_Reference_Assoc             allele/reference associations (all types)
#       MGI_Synonym                     allele synonyms
#       ALL_Allele_Mutation.bcp         Molecular mutation association
#       VOC_Annot                       allele/subtype annotations
#       All_CellLine                    mutant cell line only
#       All_Allele_CellLine             association between an Allele and a Cell Line
#
#       Diagnostics file - for verification calls to loadlib and sourceloadlib
#       Error file - for verification calls to loadlib and sourceloadlib
#
#	Annotation Load file
#
# Exit Codes:
#
# Assumes:
#
# Implementation:
#
# History
#
# 04/12/2022    sc
#       - CREAM project
#
#

import sys
import os
import db
import mgi_utils
import loadlib
import sourceloadlib

# globals

#
# from configuration file
#
inputFileName = os.getenv('INPUT_FILE_QC')
outputDir = os.getenv('OUTPUTDIR')
BCP_COMMAND = os.getenv('PG_DBUTILS') + '/bin/bcpin.csh'

# if 'true',bcp files will not be bcp-ed into the database.
# Default is 'false'
DEBUG = os.getenv('LOG_DEBUG')

#
# File descriptors
#
fpDiagFile = ''		# diagnostic file
fpErrorFile = ''	# error file 
fpInputFile = ''		

fpAlleleFile = ''       
fpMutationFile = ''	
fpRefFile = ''          
fpAccFile = ''          
fpNoteFile = ''		
fpSynonymFile = ''
fpAnnotFile = ''        
fpMutantFile = ''
fpMclFile = ''

#
# Table Names
#
alleleTable = 'ALL_Allele'
mutationTable = 'ALL_Allele_Mutation'
refTable = 'MGI_Reference_Assoc'
accTable = 'ACC_Accession'
noteTable = 'MGI_Note'
synonymTable = 'MGI_Synonym'
annotTable = 'VOC_Annot'
mclAssocTable = 'ALL_Allele_CellLine'
mclTable = 'ALL_Cellline'

#
# bcp file paths
#
alleleFileName = outputDir + '/' + alleleTable + '.bcp'
mutationFileName = outputDir + '/' + mutationTable + '.bcp'
refFileName = outputDir + '/' + refTable + '.bcp'
accFileName = outputDir + '/' + accTable + '.bcp'
noteFileName = outputDir + '/' + noteTable + '.bcp'
synonymFileName = outputDir + '/' + synonymTable + '.bcp'
annotFileName = outputDir + '/' + annotTable + '.bcp'
mclAssocFileName =  outputDir + '/' + mclAssocTable + '.bcp'
mclFileName = outputDir + '/' + mclTable + '.bcp'

#
# log file paths 
#
head, tail = os.path.split(inputFileName)

diagFileName = outputDir + '/' + tail + '.diagnostics'
errorFileName = outputDir + '/' + tail + '.error'
 
refAssocKey = 0         # MGI_Reference_Assoc._Assoc_key
accKey = 0              # ACC_Accession._Accession_key
noteKey = 0             # MGI_Note._Note_key
mgiKey = 0              # ACC_AccessionMax.maxNumericPart
annotKey = 0            # VOC_Annot._Annot_key
alleleMutationKey = 0 	# ALL_Allele_Mutation.bcp._Assoc_key

mutantAssocKey = 0	# ALL_Allele_CellLine._Assoc_key
mclKey = 0              # ALL_CellLine._CellLine_key

# Allele MGI_Note._NoteType_key 
molecularNoteTypeKey = 1021      # _NoteType_key for Molecular note
nomenNoteTypeKey = 1022          # _NoteType_key for Nomenclature note
generalNoteTypeKey = 1020        # _NoteType_key for General note
colonyIdNoteTypeKey = 1041   	 # _NoteType_key  for IKMC Allele Colony Name note

# Allele MGI_Reference_Assoc._RefAssocType_Key
origRefTypeKey = 1011           # _RefAssocType_Key for Original reference
transRefTypeKey = 1023          # _RefAssocType_Key for Transmission reference
molRefTypeKey = 1012            # _RefAssocType_Key for Molecular reference
indexRefTypeKey = 1013          # _RefAssocType_Key for Indexed reference

# Allele Synonyms
generalSynonymTypeKey = 1016    # MGI_Synonym._SynonymType_key for General synonym
synonymRefKey = 22864           # MGI_Synonym._Refs_key for J:23000 i
                                # Nomenclature Committee Use
# Allele subtypes
annotTypeKey = 1014     # Allele/Subtype
qualifierKey = 1614158  # Generic Annotation Qualifier (null term)

# Allele global Attributes
mgiTypeKey = 11
mgiPrefix = 'MGI:'
alleleKey = 0           # ALL_Allele._Allele_key
# 'Curated' from  _vocab_key = 73 (Marker-Allele Association Status)
markerStatusKey = 4268545

alleleRefKey = ''       # ALL_Allele._Refs_key - null for this load
isWildType = 0          # ALL_Allele.isWildType 
isExtinct = 0           # ALL_Allele.isExtinct
isMixed = 0             # ALL_Allele.isMixed

# MCL global attributes
esCellKey = 3982968

NS = 'Not Specified'

loaddate = loadlib.loaddate

def exit(
    status,          # numeric exit status (integer)
    message = None   # exit message (str.
    ):
    # Purpose: prints error 'message' if it is not None
    #     writes to log files and exits with 'status'
    # Returns: nothing
    # Assumes: Nothing
    # Effects: Exits with 'status'

    if message is not None:
        sys.stderr.write('\n' + str(message) + '\n')
 
    try:
        fpDiagFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        fpErrorFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        fpDiagFile.close()
        fpErrorFile.close()
        fpInputFile.close()
    except:
        pass

    db.useOneConnection(0)
    sys.exit(status)
 
def initialize():
    # Purpose: open file descriptors; write timestamps to log files
    # Returns: 1 if error, else 0
    # Assumes: Nothing
    # Effects: Sets global variables, exits if a file can't be opened,
    #  creates files in the file system

    global fpDiagFile, fpErrorFile, fpInputFile
    global fpAlleleFile, fpMutationFile, fpRefFile, fpAccFile
    global fpNoteFile, fpSynonymFile, fpAnnotFile, fpMutantFile, fpMclFile
 
    db.useOneConnection(1)
 
    try:
        fpDiagFile = open(diagFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % diagFileName)
                
    try:
        fpErrorFile = open(errorFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % errorFileName)
                
    try:
        fpInputFile = open(inputFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inputFileName)

    try:
        fpAlleleFile = open(alleleFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % alleleFileName)

    try:
        fpMutationFile = open(mutationFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % mutationFileName)

    try:
        fpRefFile = open(refFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % refFileName)

    try:
        fpAccFile = open(accFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % accFileName)

    try:
        fpNoteFile = open(noteFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % noteFileName)

    try:
        fpSynonymFile = open(synonymFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % synonymFileName)

    try:
        fpAnnotFile = open(annotFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % annotFileName)

    try:
        fpMutantFile = open(mclAssocFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % mclAssocFileName)

    try:
        fpMclFile = open(mclFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % mclFileName)

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    fpDiagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    fpDiagFile.write('Server: %s\n' % (db.get_sqlServer()))
    fpDiagFile.write('Database: %s\n' % (db.get_sqlDatabase()))

    fpErrorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

    return 0

def closeFiles():
    # Purpose: Close all file descriptors
    # Returns: 1 if error, else 0
    # Assumes: all file descriptors were initialized
    # Effects: Nothing
    # Throws: Nothing
 
    try:
        fpAlleleFile.close()
        fpMutationFile.close()
        fpRefFile.close()
        fpAccFile.close()
        fpNoteFile.close()
        fpSynonymFile.close()
        fpAnnotFile.close()
        fpMutantFile.close()
        fpMclFile.close()
    except:
        return 1
    return 0


def setPrimaryKeys():
    # Purpose: sets global primary key variables
    # Returns: 1 if error, else 0
    # Assumes: database connection exists
    # Effects: Nothing
    # Throws: Nothing

    global alleleKey, refAssocKey, accKey, noteKey, mgiKey, synonymKey, alleleMutationKey, mutantAssocKey, annotKey, mclKey

    results = db.sql(''' select nextval('all_allele_seq') as maxKey ''', 'auto')
    alleleKey = results[0]['maxKey']

    results = db.sql(''' select nextval('mgi_reference_assoc_seq') as maxKey ''', 'auto')
    refAssocKey = results[0]['maxKey']

    results = db.sql('select max(_Accession_key) + 1 as maxKey from ACC_Accession', 'auto')
    accKey = results[0]['maxKey']

    results = db.sql(''' select nextval('mgi_note_seq') as maxKey ''', 'auto')
    noteKey = results[0]['maxKey']

    results = db.sql(''' select max(maxNumericPart) + 1 as maxKey from ACC_AccessionMax where prefixPart = '%s' ''' % (mgiPrefix), 'auto')
    mgiKey = results[0]['maxKey']

    results = db.sql(''' select nextval('mgi_synonym_seq')as maxKey ''', 'auto')
    synonymKey = results[0]['maxKey']
    print('synonymKey: %s' % synonymKey)
    results = db.sql(''' select nextval('all_allele_mutation_seq') as maxKey ''', 'auto')
    alleleMutationKey = results[0]['maxKey']

    results = db.sql(''' select nextval('all_allele_cellline_seq') as maxKey ''', 'auto')
    mutantAssocKey = results[0]['maxKey']

    results = db.sql(''' select nextval('voc_annot_seq') as maxKey ''', 'auto')
    annotKey = results[0]['maxKey']

    results = db.sql(''' select nextval('all_cellline_seq') as maxKey ''', 'auto')
    mclKey = results[0]['maxKey']

    return 0

def bcpFiles():
    # Purpose: BCPs the data into the database
    # Returns: 1 if error,  else 0
    # Assumes: connection to the database
    # Effects: copies data into the db
    # Throws: Nothing

    if DEBUG  == 'true':
        return 0

    closeFiles()

    bcpI = '%s %s %s' % (BCP_COMMAND, db.get_sqlServer(), db.get_sqlDatabase())
    bcpII = '"|" "\\n" mgd'

    bcp1 = '%s %s "/" %s %s' % (bcpI, alleleTable, alleleFileName, bcpII)
    bcp2 = '%s %s "/" %s %s' % (bcpI, mutationTable, mutationFileName, bcpII)
    bcp3 = '%s %s "/" %s %s' % (bcpI, refTable, refFileName, bcpII)
    bcp4 = '%s %s "/" %s %s' % (bcpI, accTable, accFileName, bcpII)
    bcp5 = '%s %s "/" %s %s' % (bcpI, synonymTable, synonymFileName, bcpII)
    bcp6 = '%s %s "/" %s %s' % (bcpI, noteTable, noteFileName, bcpII)
    bcp7 = '%s %s "/" %s %s' % (bcpI, annotTable, annotFileName, bcpII)
    bcp8 = '%s %s "/" %s %s' % (bcpI, mclTable, mclFileName, bcpII)
    bcp9 = '%s %s "/" %s %s' % (bcpI, mclAssocTable, mclAssocFileName, bcpII)

    db.commit()

    for bcpCmd in [bcp1, bcp2, bcp3, bcp4, bcp5, bcp6, bcp7, bcp8, bcp9]:
        fpDiagFile.write('%s\n' % bcpCmd)
        os.system(bcpCmd)

    # update all_allele_seq auto-sequence
    db.sql(''' select setval('all_allele_seq', (select max(_Allele_key) from ALL_Allele)) ''', None)
    # update mgi_reference_assoc_seq auto-sequence
    db.sql(''' select setval('mgi_reference_assoc_seq', (select max(_Assoc_key) from MGI_Reference_Assoc)) ''', None)
    # update mgi_note_seq auto-sequence
    db.sql(''' select setval('mgi_note_seq', (select max(_Note_key) from MGI_Note)) ''', None)
    # update all_allele_mutation_seq auto-sequence
    db.sql(''' select setval('all_allele_mutation_seq', (select max(_Assoc_key) from ALL_Allele_Mutation)) ''', None)
    # update all_allele_cellline_seq auto-sequence
    db.sql(''' select setval('all_allele_cellline_seq', (select max(_Assoc_key) from ALL_Allele_CellLine)) ''', None)
    # update voc_annot_seq auto-sequence
    db.sql(''' select setval('voc_annot_seq', (select max(_Annot_key) from VOC_Annot)) ''', None)
    # update all_cellline_seq auto-sequence
    db.sql(''' select setval('all_cellline_seq', (select max(_CellLine_key) from ALL_CellLine)) ''', None)
    # update mgi_synonym_seq auto-sequence
    db.sql(''' select setval('mgi_synonym_seq', (select max(_Synonym_key) from MGI_Synonym)) ''', None)


    db.commit()

    return 0

def processNote(noteTypeKey, note, alleleKey, createdByKey):
    # Purpose: create note for alleleKey
    # Returns: 1 if error,  else 0
    # Assumes: file descriptor has been initialized
    # Effects: writes to the file system
    # Throws: Nothing

    global noteKey

    mgiNoteSeqNum = 1
    if note:
        fpNoteFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
            % (noteKey, alleleKey, mgiTypeKey, noteTypeKey, \
               note, createdByKey, createdByKey, loaddate, loaddate))

        noteKey += 1

    return 0

def processRefs(refTypeKey, jNums, alleleKey, createdByKey, lineNum):
    # Purpose: create reference association btwn jNum and alleleKey
    # Returns: 1 if error,  else 0
    # Assumes: file descriptor has been initialized
    # Effects: writes to the file system
    # Throws: Nothing

    global refAssocKey

    #print('jNums: %s' % jNums)
    if jNums:
        for refID in jNums.split('|'):
            refKey = loadlib.verifyReference(refID, lineNum, fpErrorFile)
            #print('refID: %s lineNum: %s refKey: %s' % (refID, lineNum, refKey))
            if refKey == 0:
                continue    # error written to fpErrorFile
            fpRefFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
                % (refAssocKey, refKey, alleleKey, mgiTypeKey, refTypeKey, \
                createdByKey, createdByKey, loaddate, loaddate))

            refAssocKey += 1
    return 0

def processSynonyms(synonyms, alleleKey, createdByKey):
    # Purpose: create synonym(s) for alleleKey
    # Returns: 1 if error,  else 0
    # Assumes: file descriptor has been initialized
    # Effects: writes to the file system
    # Throws: Nothing

    global synonymKey

    if synonyms:
        for s in synonyms.split('|'):
            fpSynonymFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
                % (synonymKey, alleleKey, mgiTypeKey, generalSynonymTypeKey, 
                synonymRefKey, s, createdByKey, createdByKey, loaddate, loaddate))

            synonymKey += 1
    return 0
    
def processSubtypes(subtypes, alleleKey, lineNum):
    # Purpose: create subtype annotations for alleleKey
    # Returns: 1 if error,  else 0
    # Assumes: file descriptor has been initialized
    # Effects: writes to the file system
    # Throws: Nothing

    global annotKey

    if subtypes:
        for s in subtypes.split('|'):
            #print('lineNum: %s subtype: %s' % (lineNum, s))
            # _vocab_key = 93 (Allele Subtype)
            alleleSubtypeKey = loadlib.verifyTerm('', 93, s, lineNum, fpErrorFile)

            fpAnnotFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
                % (annotKey, annotTypeKey, alleleKey, alleleSubtypeKey, \
                    qualifierKey, loaddate, loaddate))
            annotKey += 1
    return 0

def processMutations(molMuts, alleleKey, lineNum):
    global alleleMutationKey

    for m in molMuts.split('|'):
            #print('lineNum: %s mutation: %s' % (lineNum, m))
            # _vocab_key = 36 (Allele Molecular Mutation)
            mutationTermKey = loadlib.verifyTerm('', 36, m, lineNum, fpErrorFile)
            fpMutationFile.write('%s|%s|%s|%s|%s\n' \
                % (alleleMutationKey, alleleKey, mutationTermKey, loaddate, loaddate))
            alleleMutationKey += 1
    return 0

def processMCLs(mclKeyList, derivationKey, strainOfOriginKey, alleleKey, createdByKey):
    # Purpose: create NS MCL and/or MCL association to allele. 
    # Returns: 1 if error,  else 0
    # Assumes: file descriptors have been initialized. mclKeyList or derivationKeyList is not empty
    # Effects: writed to the file system
    # Throws: Nothing
    global mutantAssocKey, mclKey

    if mclKeyList:
        # create MCL association(s) to the allele
        for m in str.split(mclKeyList, '|'):
            fpMutantFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
                % (mutantAssocKey, alleleKey, m, \
                    createdByKey, createdByKey, loaddate, loaddate))
            mutantAssocKey += 1

    else:
        # otherwise create new not specified MCL with derivation
        # and MCL association to the allele
        fpMclFile.write('%s|%s|%s|%s|%s|1|%s|%s|%s|%s\n' \
            % (mclKey, NS, esCellKey, strainOfOriginKey, derivationKey, \
                createdByKey, createdByKey, loaddate, loaddate))

        fpMutantFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
            % (mutantAssocKey, alleleKey, mclKey, \
                createdByKey, createdByKey, loaddate, loaddate))
        mclKey += 1
        mutantAssocKey += 1

    return 0
def processFile():
    # Purpose: Read the input file, resolve values to keys. Create bcp files
    # Returns: 1 if error,  else 0
    # Assumes: file descriptors have been initialized
    # Effects: exits if the line does not have 23 columns
    # Throws: Nothing

    global alleleKey, accKey, mgiKey 

    lineNum = 0
    # For each line in the input file

    for line in fpInputFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = line[:-1].split('\t')
        
        try:
            aSymbol = tokens[0]
            aName= tokens[1]
            geneID = tokens[2]
            user = tokens[3]
            aStatus = tokens[4]
            aType = tokens[5]
            inheritMode = tokens[6]
            transmission = tokens[7]
            collection = tokens[8]
            molNote = tokens[9]
            nomenNote = tokens[10]
            genNotes = tokens[11]
            colonyNote = tokens[12]
            origRef = tokens[13]
            transRef = tokens[14]
            molRef = tokens[15]
            idxRefs = tokens[16]
            synonyms = tokens[17]
            subtypes = tokens[18]
            molMuts = tokens[19]
            pcl = tokens[20]
            soo = tokens[21]
            mclKeyList = tokens[22]
            derivationKey = tokens[23]
            
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

        # marker key
        markerKey = loadlib.verifyMarker(geneID, lineNum, fpErrorFile)
            
        # creator
        createdByKey = loadlib.verifyUser(user, lineNum, fpErrorFile)

        # _vocab_key = 37 (Allele Status)
        alleleStatusKey = loadlib.verifyTerm('', 37, aStatus, lineNum, fpErrorFile)

        # _vocab_key = 38 (Allele Type)
        alleleTypeKey = loadlib.verifyTerm('', 38, aType, lineNum, fpErrorFile)

        # _vocab_key = 35 (Allele Inheritance Mode)
        inheritanceModeKey = loadlib.verifyTerm('', 35, inheritMode, lineNum, fpErrorFile)

        # _vocab_key = 61 (Allele Transmission)
        transmissionKey = loadlib.verifyTerm('', 61, transmission, lineNum, fpErrorFile)

        # _vocab_key = 92 (Allele Collection)
        collectionKey = loadlib.verifyTerm('', 92, collection, lineNum, fpErrorFile)

        # strain of origin
        strainOfOriginKey = sourceloadlib.verifyStrain(soo, lineNum, fpErrorFile)

        # if errors, continue to next record
        # errors are stored (via loadlib) in the .error log

        if markerKey == 0 \
                or createdByKey == 0 \
                or alleleStatusKey == 0 \
                or alleleTypeKey == 0 \
                or inheritanceModeKey == 0 \
                or transmissionKey == 0 \
                or collectionKey == 0 \
                or strainOfOriginKey == 0:
            continue

        # if no errors, process the allele

        # allele (master)
        fpAlleleFile.write('%d|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
            % (alleleKey, markerKey, strainOfOriginKey, inheritanceModeKey,  \
            alleleTypeKey, alleleStatusKey, transmissionKey, collectionKey, aSymbol,\
            aName, isWildType, isExtinct, isMixed, alleleRefKey, markerStatusKey, \
            createdByKey, createdByKey, createdByKey, loaddate, loaddate, \
            loaddate))

        # MGI ID for the llele
        fpAccFile.write('%s|%s%d|%s|%s|1|%d|%d|0|1|%s|%s|%s|%s\n' \
            % (accKey, mgiPrefix, mgiKey, mgiPrefix, mgiKey, alleleKey, mgiTypeKey, \
               createdByKey, createdByKey, loaddate, loaddate))

        # process the Notes
        processNote(molecularNoteTypeKey, molNote, alleleKey, createdByKey)
        processNote(nomenNoteTypeKey, nomenNote, alleleKey, createdByKey)
        processNote(generalNoteTypeKey, genNotes, alleleKey, createdByKey)
        processNote(colonyIdNoteTypeKey, colonyNote, alleleKey, createdByKey)

        # process the references
        #print('original ref')
        processRefs(origRefTypeKey, origRef, alleleKey, createdByKey, lineNum)
        #print('transmissiont ref')
        processRefs(transRefTypeKey, transRef, alleleKey, createdByKey, lineNum)
        #print('molecular ref')
        processRefs(molRefTypeKey, molRef, alleleKey, createdByKey, lineNum)
        #print('indexed ref')
        processRefs(indexRefTypeKey, idxRefs, alleleKey, createdByKey, lineNum)

        # process synonyms, subtypes, mutations and mcls
        processSynonyms(synonyms, alleleKey, createdByKey)
        processSubtypes(subtypes, alleleKey, lineNum)
        processMutations(molMuts, alleleKey, lineNum)

        # if either mclKeyList or derivationKey have data, then process the MCLs
        if mclKeyList or derivationKey:
            print('mclKeyList: %s derivationKey: %s' % (mclKeyList, derivationKey))
            processMCLs(mclKeyList, derivationKey, strainOfOriginKey, alleleKey, createdByKey)
 
        accKey += 1
        mgiKey += 1
        alleleKey += 1

    #   end of "for line in inputFile.readlines():"

    #
    # Update the AccessionMax value
    #

    if not DEBUG:
        db.sql('select * from ACC_setMax(%d)' % (lineNum), None)
        db.commit()

    return 0

#
# MAIN
#
rc = 0
if initialize() != 0:
    rc = 1

if setPrimaryKeys() != 0:
    rc = 1

if processFile() != 0:
    rc = 1

if bcpFiles() != 0:
    rc = 1

fpDiagFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
fpDiagFile.close()
fpErrorFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
fpErrorFile.close()

sys.exit(rc)

