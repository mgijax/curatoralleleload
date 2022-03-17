#
# makeAllele.py
###############################################################################
#
# Purpose:
#
#	A generic allele load to load alleles and their notes, references, 
#           synonyms, mutations and subtypes
#
# Usage:
#	makeAllele.py 
#
# Envvars:
#
# Inputs:
#
#	A tab-delimited file in the format:
#
#	field 1:  MGI Marker ID
#	field 2:  Strain of Origin Name
#       field 3:  Inheritance Mode
#	field 4:  Allele Type
#	field 5:  Allele Status
#	field 6:  Transmission
#       field 7:  Allele Collection
#       field 8:  Allele Symbol
#       field 9:  Allele Name
#       field 10: Created By
#       field 11: Notes - all types, multivalued 'type1||note|||type2||note'
#       field 12: References - JNum all types, multivalued 'type1||jnum|||type2||jnum'
#       field 13: Synonyms - all types, multivalued 'type1||syn1|||type2||syn2'
#       field 14: Molecular Mutation, multivalued '||' delimitted
#       field 15: Allele Subtype, multivalued '||' delimitted
#
# Outputs:
#
#       BCP files:
#
#       ALL_Allele.bcp                  master Allele records
#       ACC_Accession.bcp               Accession records
#       ACC_AccessionReference.bcp      Accession Reference records
#       MGI_Note                        allele notes (all types)
#       MGI_Reference_Assoc             allele/reference associations (all types)
#       MGI_Synonym                     allele synonyms
#       ALL_Allele_Mutation.bcp         Molecular mutation association
#       VOC_Annot                       allele/subtype annotations
#
#       Diagnostics file of all input parameters and SQL commands
#       Error file
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
# 02/22/2022    sc
#       - CREAM project
#
#

import sys
import os
import db
import mgi_utils
import loadlib
import sourceloadlib

#globals

#
# from configuration file
#
user = os.environ['MGD_DBUSER']
passwordFileName = os.environ['MGD_DBPASSWORDFILE']
inputFileName = os.environ['INPUTFILE']
outputDir = os.environ['OUTPUTDIR']
jnum = os.environ['JNUMBER']
BCP_COMMAND = os.environ['PG_DBUTILS'] + '/bin/bcpin.csh'

DEBUG = 0		# if 0, not in debug mode

bcpon = 1		# can the bcp files be bcp-ed into the database?  default is yes.

diagFile = ''		# diagnostic file descriptor
errorFile = ''		# error file descriptor
inputFile = ''		# file descriptor
alleleFile = ''         # file descriptor
mutationFile = ''	# file descriptor
mutantFile = ''		# file descriptor
refFile = ''            # file descriptor
accFile = ''            # file descriptor
accRefFile = ''         # file descriptor
noteFile = ''		# file descriptor
annotFile = ''		# file descriptor
newAlleleFile = ''      # file descriptor

alleleTable = 'ALL_Allele'
mutationTable = 'ALL_Allele_Mutation'
mutantTable = 'ALL_Allele_CellLine'
refTable = 'MGI_Reference_Assoc'
accTable = 'ACC_Accession'
accRefTable = 'ACC_AccessionReference'
noteTable = 'MGI_Note'
annotTable = 'VOC_Annot'

alleleFileName = outputDir + '/' + alleleTable + '.bcp'
mutationFileName = outputDir + '/' + mutationTable + '.bcp'
mutantFileName =  outputDir + '/' + mutantTable + '.bcp'
refFileName = outputDir + '/' + refTable + '.bcp'
accFileName = outputDir + '/' + accTable + '.bcp'
accRefFileName = outputDir + '/' + accRefTable + '.bcp'
noteFileName = outputDir + '/' + noteTable + '.bcp'
annotFileName = outputDir + '/' + annotTable + '.bcp'

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name
newAlleleFileName = ''	# output file with new accession ids

alleleKey = 0           # ALL_Allele._Allele_key
mutantionKey = 0 	# ALL_Allele_Mutation.bcp._Assoc_key
mutantKey = 0  		# ALL_Allele_CellLine._Assoc_key
refAssocKey = 0		# MGI_Reference_Assoc._Assoc_key
accKey = 0              # ACC_Accession._Accession_key
noteKey = 0		# MGI_Note._Note_key
mgiKey = 0              # ACC_AccessionMax.maxNumericPart
annotKey = 0		# VOC_Annot._Annot_key
mgiNoteObjectKey = 11   # MGI_Note._MGIType_key
mgiMolecularNoteTypeKey = 1021   # MGI_Note._NoteType_key
mgiDriverNoteTypeKey = 1034   	 # MGI_Note._NoteType_key
mgiIKMCNoteTypeKey = 1041   	 # MGI_Note._NoteType_key
ikmcSQLs = []

mgiTypeKey = 11		# Allele
mgiPrefix = 'MGI:'
annotTypeKey = 1014
qualifierKey = 1614158

# key = symbol
# value = (alleleKey, noteKey, mgiKey)
alleleLookup = {}

loaddate = loadlib.loaddate

#
# Purpose: prints error message and exits
#
def exit(
    status,          # numeric exit status (integer)
    message = None   # exit message (str.
    ):

    if message is not None:
        sys.stderr.write('\n' + str(message) + '\n')
 
    try:
        diagFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        errorFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        diagFile.close()
        errorFile.close()
        inputFile.close()
        newAlleleFile.close()
    except:
        pass

    db.useOneConnection(0)
    sys.exit(status)
 
#
# Purpose: process command line options
#
def initialize():
    global diagFile, errorFile, inputFile, errorFileName, diagFileName
    global alleleFile, mutationFile, mutantFile, refFile
    global accFile, accRefFile, noteFile, annotFile
    global newAlleleFile
 
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)
 
    head, tail = os.path.split(inputFileName) 

    diagFileName = outputDir + '/' + tail + '.diagnostics'
    errorFileName = outputDir + '/' + tail + '.error'
    newAlleleFileName = outputDir + '/' + tail + '.new'

    try:
        diagFile = open(diagFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % diagFileName)
                
    try:
        errorFile = open(errorFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % errorFileName)
                
    try:
        newAlleleFile = open(newAlleleFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % newAlleleFileName)

    try:
        inputFile = open(inputFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inputFileName)

    try:
        alleleFile = open(alleleFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % alleleFileName)

    try:
        mutationFile = open(mutationFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % mutationFileName)

    try:
        mutantFile = open(mutantFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % mutantFileName)

    try:
        refFile = open(refFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % refFileName)

    try:
        accFile = open(accFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % accFileName)

    try:
        accRefFile = open(accRefFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % accRefFileName)

    try:
        noteFile = open(noteFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % noteFileName)

    try:
        annotFile = open(annotFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % annotFileName)

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.write('Server: %s\n' % (db.get_sqlServer()))
    diagFile.write('Database: %s\n' % (db.get_sqlDatabase()))

    errorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

#
# Purpose: Close files.
#
def closeFiles():

    alleleFile.close()
    mutationFile.close()
    mutantFile.close()
    refFile.close()
    accFile.close()
    accRefFile.close()
    noteFile.close()
    annotFile.close()

#
# Purpose:  sets global primary key variables
#
def setPrimaryKeys():

    global alleleKey, refAssocKey, accKey, noteKey, mgiKey, mutationKey, mutantKey, annotKey

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

    results = db.sql(''' select nextval('all_allele_mutation_seq') as maxKey ''', 'auto')
    mutationKey = results[0]['maxKey']

    results = db.sql(''' select nextval('all_allele_cellline_seq') as maxKey ''', 'auto')
    mutantKey = results[0]['maxKey']

    results = db.sql(''' select nextval('voc_annot_seq') as maxKey ''', 'auto')
    annotKey = results[0]['maxKey']

#
# Purpose:  BCPs the data into the database
#
def bcpFiles():

    bcpdelim = "|"

    if DEBUG or not bcpon:
        return

    closeFiles()

    bcpI = '%s %s %s' % (BCP_COMMAND, db.get_sqlServer(), db.get_sqlDatabase())
    bcpII = '"|" "\\n" mgd'

    bcp1 = '%s %s "/" %s %s' % (bcpI, alleleTable, alleleFileName, bcpII)
    bcp2 = '%s %s "/" %s %s' % (bcpI, mutationTable, mutationFileName, bcpII)
    bcp3 = '%s %s "/" %s %s' % (bcpI, mutantTable, mutantFileName, bcpII)
    bcp4 = '%s %s "/" %s %s' % (bcpI, refTable, refFileName, bcpII)
    bcp5 = '%s %s "/" %s %s' % (bcpI, accTable, accFileName, bcpII)
    bcp6 = '%s %s "/" %s %s' % (bcpI, accRefTable, accRefFileName, bcpII)
    bcp7 = '%s %s "/" %s %s' % (bcpI, noteTable, noteFileName, bcpII)
    bcp8 = '%s %s "/" %s %s' % (bcpI, annotTable, annotFileName, bcpII)

    db.commit()

    for bcpCmd in [bcp1, bcp2, bcp3, bcp4, bcp5, bcp6, bcp7, bcp8]:
        diagFile.write('%s\n' % bcpCmd)
        os.system(bcpCmd)

    if len(ikmcSQLs) > 0:
        print(ikmcSQLs)
        db.sql(ikmcSQLs, None)
        db.commit()

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

    db.commit()

#
# Purpose:  processes data
#
# a) add additional mutant cell lines to new and existing alleles
# b) add additional IKMC Colony note to new and existing alleles
# c) set Allele Status = Approved for reserved alleles
#
def processFileIKMC(createMCL, createNote, setStatus, \
        symbol, ikmcSymbol, mutantCellLine, ikmcNotes, createdByKey, existingAlleleID):

    global noteKey, ikmcSQLs

    #
    # add new MCLs to new/existing alleles
    #
    if len(createMCL) > 0:

        if DEBUG:
                print(symbol, createMCL)

        if int(createMCL) == 0:
                aKey = alleleLookup[symbol][0][0]
        else:
                aKey = createMCL

        addMutantCellLine(aKey, mutantCellLine, createdByKey)

    #
    # set allele/status = Approved for existing "reserved" alleles
    #
    if len(setStatus) > 0:
        ikmcSQLs.append('update ALL_Allele set _Allele_Status_key = 847114 where _Allele_key = %s' % (setStatus))

    #
    # Add IKMC Colony/Note to a new or existing allele
    #
    # child exists/ikmc note exists : update existing note
    # 	|| => _Note_key||existing colony notes
    #
    # child exists/ikmc note does not exis : add note
    # 	:: => allele/child key
    #
    # new allele/child/non-duplicate IKMC Colony
    #	0::colony(s)
    #
    # blank => do nothing
    #

    if len(createNote) > 0:

        if DEBUG:
                print('createNote: ', symbol)

        try:
            tokens = createNote.split('::')
            aKey = tokens[0]

            # duplicate child, additional note : add note to new child
            if int(aKey) == 0:
                nKey = alleleLookup[symbol][0][1]
                note = tokens[1]
                ikmcSQLs.append('''update MGI_Note set note = '%s' where _Note_key = %s;''' % (note, nKey))
                        
            # child exists, note does not exist : add note to existing child
            else:
                aKey = tokens[0]
                note = ikmcNotes

                if symbol in alleleLookup:
                        nKey = alleleLookup[symbol][0][1]
                        ikmcSQLs.append('''update MGI_Note set note = rtrim(note) || '|%s' where _Note_key = %s;''' % (note, nKey))
                else:
                        noteFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
                        % (noteKey, aKey, mgiNoteObjectKey, mgiIKMCNoteTypeKey, \
                        note, createdByKey, createdByKey, loaddate, loaddate))

                        # save symbol/aKey/ikmc note key/allele id
                        alleleLookup[symbol] = []
                        alleleLookup[symbol].append((aKey, noteKey, 'missing allele id (1)'))

                        noteKey = noteKey + 1

        # child exists, note exists : update existing note
        except:
            if DEBUG:
                print(createNote)

            tokens = createNote.split('||')
            nKey = tokens[0]
            note = tokens[1] + '|' + ikmcNotes
            ikmcSQLs.append('''update MGI_Note set note = '%s' where _Note_key = %s;''' % (note, nKey))
                        
    # 
    # print out the proper allele id
    #
    if len(existingAlleleID) > 0:
        printAlleleID = existingAlleleID
    elif symbol in alleleLookup:
        printAlleleID = alleleLookup[symbol][0][2]
    else:
        printAlleleID = 'missing allele id (2)'

    newAlleleFile.write('%s\t%s\t%s\n' \
                % (mgi_utils.prvalue(ikmcNotes), \
                        mgi_utils.prvalue(printAlleleID), \
                        mgi_utils.prvalue(ikmcSymbol)))

#
# Purpose:  processes data
#
def processFile():

    global alleleKey, refAssocKey, accKey, noteKey, mgiKey, annotKey, mutationKey
    global alleleLookup

    lineNum = 0
    # For each line in the input file

    for line in inputFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = line[:-1].split('\t')
        #print line
        try:
            markerID = tokens[0]
            symbol = tokens[1]
            name = tokens[2]
            alleleStatus = tokens[3]
            alleleType = tokens[4]
            alleleSubtypes = tokens[5]
            collectionKey = tokens[6]
            germLine = tokens[7]
            references = tokens[8]
            strainOfOrigin = tokens[9]
            mutantCellLine = tokens[10]
            molecularNotes = tokens[11]
            driverNotes = tokens[12]
            ikmcNotes = tokens[13]
            mutations = tokens[14]
            inheritanceMode = tokens[15]
            isMixed = tokens[16]
            isExtinct = tokens[17]
            createdBy = tokens[18]
            createMCL = tokens[19]
            createNote = tokens[20]
            setStatus = tokens[21]
            existingAlleleID = tokens[22]
            ikmcSymbol = tokens[23]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

        # creator
        createdByKey = loadlib.verifyUser(createdBy, lineNum, errorFile)
        if createdByKey == 0:
            continue

        # processing for IKMC-only
        if len(createMCL) > 0 or len(createNote) > 0 or len(setStatus) > 0:
                processFileIKMC(createMCL, createNote, setStatus, \
                        symbol, ikmcSymbol, mutantCellLine, ikmcNotes, \
                        createdByKey, existingAlleleID)
                continue

        # marker key
        markerKey = loadlib.verifyMarker(markerID, lineNum, errorFile)

        # hard-coded
        # _vocab_key = 73 (Marker-Allele Association Status)
        # _term_key = 4268545 (Curated)
        markerStatusKey = 4268545

        # _vocab_key = 37 (Allele Status)
        alleleStatusKey = loadlib.verifyTerm('', 37, alleleStatus, lineNum, errorFile)

        # _vocab_key = 38 (Allele Type)
        alleleTypeKey = loadlib.verifyTerm('', 38, alleleType, lineNum, errorFile)

        # _vocab_key = 61 (Allele Transmission)
        germLineKey = loadlib.verifyTerm('', 61, germLine, lineNum, errorFile)

        # _vocab_key = 36 (Allele Molecular Mutation)
        allMutations = mutations.split('|')

        # _vocab_key = 35 (Allele Status)
        inheritanceModeKey = loadlib.verifyTerm('', 35, inheritanceMode, lineNum, errorFile)

        # strains
        strainOfOriginKey = sourceloadlib.verifyStrain(strainOfOrigin, lineNum, errorFile)

        # reference
        refKey = loadlib.verifyReference(jnum, lineNum, errorFile)

        # if errors, continue to next record
        # errors are stored (via loadlib) in the .error log

        if markerKey == 0 \
                or markerStatusKey == 0 \
                or alleleStatusKey == 0 \
                or alleleTypeKey == 0 \
                or germLineKey == 0 \
                or allMutations == 0 \
                or inheritanceModeKey == 0 \
                or strainOfOriginKey == 0 \
                or refKey == 0 \
                or createdByKey == 0:
            continue

        # if no errors, process the allele

        # not specified/testing
        #collectionKey = 11025586

        # allele (master)
        alleleFile.write('%d|%s|%s|%s|%s|%s|%s|%s|%s|%s|0|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
            % (alleleKey, markerKey, strainOfOriginKey, inheritanceModeKey, alleleTypeKey, \
            alleleStatusKey, germLineKey, collectionKey, symbol, name, \
            isExtinct, isMixed, refKey, markerStatusKey, \
            createdByKey, createdByKey, createdByKey, loaddate, loaddate, loaddate))

        # molecular mutation
        for mutation in allMutations:
                mutationTermKey = loadlib.verifyTerm('', 36, mutation, lineNum, errorFile)
                mutationFile.write('%s|%s|%s|%s|%s\n' \
                % (mutationKey, alleleKey, mutationTermKey, loaddate, loaddate))
                mutationKey = mutationKey + 1

        #
        # allele references
        #
        allReferences = references.split('||')
        for reference in allReferences:
                refType, refID = reference.split('|')
                refKey = loadlib.verifyReference(refID, lineNum, errorFile)

                if refType == 'Original':
                        refAssocTypeKey = 1011
                elif refType == 'Transmission':
                        refAssocTypeKey = 1023
                elif refType == 'Molecular':
                        refAssocTypeKey = 1012

                refFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
                        % (refAssocKey, refKey, alleleKey, mgiTypeKey, refAssocTypeKey, \
                        createdByKey, createdByKey, loaddate, loaddate))
                refAssocKey = refAssocKey + 1

        #
        # allele subtypes
        #
        allSubtypes = alleleSubtypes.split('|')
        for s in allSubtypes:

                # _vocab_key = 93 (Allele Subtype)
                alleleSubtypeKey = loadlib.verifyTerm('', 93, s, lineNum, errorFile)

                annotFile.write('%s|%s|%s|%s|%s|%s|%s\n' \
                        % (annotKey, annotTypeKey, alleleKey, alleleSubtypeKey, \
                                qualifierKey, loaddate, loaddate))
                annotKey = annotKey + 1

        #
        # mutant cell line
        #
        if len(mutantCellLine) > 0:
            addMutantCellLine(alleleKey, mutantCellLine, createdByKey)

        # MGI Accession ID for the allelearker

        accFile.write('%s|%s%d|%s|%s|1|%d|%d|0|1|%s|%s|%s|%s\n' \
            % (accKey, mgiPrefix, mgiKey, mgiPrefix, mgiKey, alleleKey, mgiTypeKey, \
               createdByKey, createdByKey, loaddate, loaddate))

        # storing data in MGI_Note
        # molecular notes

        mgiNoteSeqNum = 1
        if len(molecularNotes) > 0:

            noteFile.write('%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
                % (noteKey, alleleKey, mgiNoteObjectKey, mgiMolecularNoteTypeKey, \
                   molecularNotes, createdByKey, createdByKey, loaddate, loaddate))

            noteKey = noteKey + 1

        # driver notes
        # TR12662/MGI_Relationship._Category_key = 1006
        # removed noteFile code
        # place hodler for MGI_Relationship code
        # the IKMC is the only product using this and IKMC does not add any driver note
        #mgiNoteSeqNum = 1
        #if len(driverNotes) > 0:

        # ikmc notes
        useIKMCnotekey = 0
        if len(ikmcNotes) > 0:

            noteFile.write('%s|%s|%s|%s|%s|%s|%s|%s\n' \
                % (noteKey, alleleKey, mgiNoteObjectKey, mgiIKMCNoteTypeKey, \
                   ikmcNotes, createdByKey, createdByKey, loaddate, loaddate))

            useIKMCnotekey = noteKey
            noteKey = noteKey + 1

        # Print out a new text file and attach the new MGI Allele IDs as the last field

        if createdBy == 'ikmc_alleleload':
                newAlleleFile.write('%s\t%s%s\t%s\n' \
                % (mgi_utils.prvalue(ikmcNotes), \
                        mgi_utils.prvalue(mgiPrefix), mgi_utils.prvalue(mgiKey), \
                        mgi_utils.prvalue(ikmcSymbol)))
        else:
                newAlleleFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s\n' \
                % (mgi_utils.prvalue(markerID), \
                mgi_utils.prvalue(symbol), \
                mgi_utils.prvalue(name), \
                mgi_utils.prvalue(alleleStatus), \
                mgi_utils.prvalue(alleleType), \
                mgi_utils.prvalue(alleleSubtype), \
                mgi_utils.prvalue(collection), \
                mgi_utils.prvalue(germLine), \
                mgi_utils.prvalue(references), \
                mgi_utils.prvalue(strainOfOrigin), \
                mgi_utils.prvalue(mutantCellLine), \
                mgi_utils.prvalue(allMutations), \
                mgi_utils.prvalue(inheritanceMode), \
                mgi_utils.prvalue(isMixed), \
                mgi_utils.prvalue(isExtinct), \
                mgi_utils.prvalue(refKey), \
                mgi_utils.prvalue(markerStatusKey), \
                mgi_utils.prvalue(createdBy), \
                mgi_utils.prvalue(mgiPrefix), mgi_utils.prvalue(mgiKey)))

        # save symbol/alleleKey/ikmc note key
        alleleLookup[symbol] = []
        alleleLookup[symbol].append((alleleKey, useIKMCnotekey, mgiPrefix + str(mgiKey)))

        accKey = accKey + 1
        mgiKey = mgiKey + 1
        alleleKey = alleleKey + 1

    #	end of "for line in inputFile.readlines():"

    #
    # Update the AccessionMax value
    #

    if not DEBUG:
        db.sql('select * from ACC_setMax(%d)' % (lineNum), None)
        db.commit()


#
# Purpose: write 1 or more mutation cell line associations to bcp file
#
def addMutantCellLine(alleleKey, mutantCellLine, createdByKey):

    global mutantKey

    mutantCellLineKey = 0

    results = db.sql('''
        select _CellLine_key from ALL_CellLine
        where isMutant = 1 and _Derivation_key is not null
        and cellLine = '%s'
        ''' % (mutantCellLine) , 'auto')

    for r in results:
        mutantCellLineKey = r['_CellLine_key']

    mutantFile.write('%d|%s|%s|%s|%s|%s|%s\n' \
                % (mutantKey, alleleKey, mutantCellLineKey, \
                createdByKey, createdByKey, loaddate, loaddate))

    mutantKey = mutantKey + 1

#
# Main
#

if __name__ == '__main__':

        initialize()

        setPrimaryKeys()

        processFile()

        bcpFiles()
