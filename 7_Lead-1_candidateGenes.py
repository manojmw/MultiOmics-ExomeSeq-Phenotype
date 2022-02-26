#!/usr/bin/python

import argparse, sys
import pandas as pd
import logging
import time

###########################################################

# Parses tab-seperated canonical transcripts file
# Required columns are: 'ENSG' and 'GENE' (can be in any order,
# but they MUST exist)
#
# Returns a dictionary:
# Key -> Gene
# Value -> ENSG
def ENSG_Gene(inCanonicalFile):

    ENSG_Gene_dict = {} # Initializing an empty dictionary

    Canonical_File = open(inCanonicalFile)

    Canonical_header_line = Canonical_File.readline() # Grabbing the header line

    Canonical_header_fields = Canonical_header_line.split('\t')

    # Check the column headers and grab indexes of our columns of interest
    (ENSG_col, Gene_col) = (-1,-1)

    for i in range(len(Canonical_header_fields)):
        if Canonical_header_fields[i] == 'ENSG':
            ENSG_col  = i
        elif Canonical_header_fields[i] == 'GENE':
            Gene_col = i

    if not ENSG_col >= 0:
        sys.exit("Missing required column title 'ENSG' in the file: %s \n" % inCanonicalFile)
    elif not Gene_col >= 0:
        sys.exit("Missing required column title 'GENE' in the file: %s \n" % inCanonicalFile)
    # else grabbed the required column indexes -> PROCEED

    # Parsing the Uniprot Primary Accession file
    for line in Canonical_File:
        line = line.rstrip('\n')
        CanonicalTranscripts_fields = line.split('\t')

        # Key ->  Gene
        # Value -> ENSG
        ENSG_Gene_dict[CanonicalTranscripts_fields[Gene_col]] = CanonicalTranscripts_fields[ENSG_col]

    # Closing the file
    Canonical_File.close()

    return ENSG_Gene_dict

###########################################################

# Parses the candidateGenes file in .xlsx format
# Required columns are: 'Gene', 'pathologyID' and
# 'Confidence score' (can be in any order,
# but they MUST exist)
#
# Returns a list with 3 items (in each sublist:
# - Gene
# - pathologyID
# - Confidence score)
def CandidateGeneParser(inCandidateFile):

    # Input - List of candidate gene file(s)
    candidate_files = inCandidateFile

    # Initializing an empty data frame to store the data from all files
    meta_data = pd.DataFrame()

    # Iterating over the list of files and appending to the DataFrame (meta_data)
    for file in candidate_files:
        data = pd.read_excel(file)
        meta_data = pd.concat([meta_data, data])

    # Extract Gene, pathologyID and Confidence score and drop rows with missing values
    CandidateGene_data = pd.DataFrame(meta_data, columns=['Gene', 'pathologyID', 'Confidence score']).dropna()

    return CandidateGene_data.values.tolist()

###########################################################

# Parses the dictionary (ENSG_Gene_dict) returned
# by the function: ENSG_Gene and the list (CandidateGene_data)
# returned by the function: CandidateGeneParser
# Calls the functions ENSG_Gene and CandidateGeneParser
#
# Maps Gene name to ENSG
# Returns a list with 3 items (in each sublist:
# - ENSG
# - pathologyID
# - Confidence score)
def CandidateGene2ENSG(inCanonicalFile, inCandidateFile):

    # Calling the function ENSG_Gene
    ENSG_Gene_dict = ENSG_Gene(inCanonicalFile)

    # Calling the function CandidateGeneParser
    CandidateGene_data = CandidateGeneParser(inCandidateFile)

    # Initializing output list
    canidateGene_out_list = []

    # Counter for canidate genes not found in the canonical transcripts file
    lost_CandidateGene = 0

    for data in CandidateGene_data:
        if data[0] in ENSG_Gene_dict.keys():
            canidateGene_out = [ENSG_Gene_dict.get(data[0]),  data[1],  data[2]]
            canidateGene_out_list.append(canidateGene_out)
        else:
            lost_CandidateGene += 1

    logging.debug("\nNo. of candidate genes not found in the Canonical transcripts file: %d" % lost_CandidateGene)

    return canidateGene_out_list

###########################################################

# Parses the High-quality Interactome produced by Interactome.py
# Required columns:
# First column - ENSG of Protein A
# Second column - ENSG of Protein B
#
# Returns a list with 2 items (in each sub list):
# - ENSG of Protein A
# - ENSG of Protein B
def Interacting_Proteins(inInteractome):

    # Input - Interactome file
    Interactome_File = open(inInteractome)

    # Initializing Interactome dictionary
    Interactome_list = []

    for line in Interactome_File:
        line = line.rstrip('\n')

        Interactome_fields = line.split('\t')

        Interacting_Proteins = [Interactome_fields[0], Interactome_fields[1]]

        Interactome_list.append(Interacting_Proteins)

    # Closing the file
    Interactome_File.close()

    return Interactome_list

###########################################################

# Parses the canidateGene_out_list returned by the function: CandidateGene2ENSG
# Checks the Number of interactors for each candidate gene
# using the Interactome_list returned by the function: Interacting_Proteins
# Checks the number of Interactors that are known candidate genes
# Returns a list

def Lead1_CandidateENSG(inCanonicalFile, inCandidateFile, inInteractome):

    # Calling the functions
    canidateGene_out_list = CandidateGene2ENSG(inCanonicalFile, inCandidateFile)
    Interactome_list = Interacting_Proteins(inInteractome)

    # Dictionary for storing the Candidate Genes and Interactors
    candGene_Interactors_list = []

    # Keep the count of candidate genes not interacting with any protein
    nonInteracting_candGene_count = 0

    # List of non-interacting candidate genes
    nonInteracting_candGeneList = []

    for candidateGene in canidateGene_out_list:

        # List for interacting proteins
        Interactors = []

        # List for interacting proteins that are known candidate genes
        Known_interactor = []

        for Proteins in Interactome_list:
            # If candidate gene is protein A
            if (candidateGene[0] == Proteins[0]):
                # then, get the interacting protein (protein B)
                Interactors.append(Proteins[1])

            # If candidate gene is protein B
            elif (candidateGene[0] == Proteins[1]):
                # then, get the interacting protein (protein A)
                Interactors.append(Proteins[0])

        candGene_Interactors = [candidateGene[0], candidateGene[1], candidateGene[2], str(len(Interactors)), Interactors]
        candGene_Interactors_list.append(candGene_Interactors)

    # Checking the number of interactors that are known canidate genes
    for candidateGene in canidateGene_out_list:
        for data in candGene_Interactors_list:
            for interactor in data[4]:
                if interactor in candidateGene:
                    Known_interactor.append(interactor)
                    data.append(Known_interactor)

    return candGene_Interactors_list

###########################################################

# Parses the list (candGene_Interactors_list)
# returned by the function: Lead1_CandidateENSG
# Maps the ENSGs in the list to Gene names using the dictionary: ENSG_Gene_dict
#
# Prints to STDOUT in .tsv format
# Output consists of 6 columns:
# - Candidate Gene
# - Pathology Identifier
# - Confidence Score
# - Count of Proteins interacting with Candidate Gene
# - Count of Interacting genes that are known candidate Genes
# - List of Known Interacting genes
def Lead1_CandidateGene(args):

    candGene_Interactors_list = Lead1_CandidateENSG(args.inCanonicalFile, args.inCandidateFile, args.inInteractome)
    ENSG_Gene_dict = ENSG_Gene(args.inCanonicalFile)

    # Printing the header for the output
    header = ['Candidate_Gene', 'pathologyID', 'Confidence_Score', 'No_of_Interactors', 'No_of_KnownInteractors', 'List_KnownInteractors']
    print('\t'.join(header))

    for CandidateENSG_data in candGene_Interactors_list:

        # known interactor genes
        Known_Interactors_Genes = []

        candidateGene = [GeneName for (GeneName, ENSG) in ENSG_Gene_dict.items() if ENSG == CandidateENSG_data[0]]
        if len(CandidateENSG_data) > 5:
            for Interactors_ENSG in CandidateENSG_data[5]:
                Known_Interactor_GeneList = [GeneName for (GeneName, ENSG) in ENSG_Gene_dict.items() if ENSG == Interactors_ENSG]
                Known_Interactor_GeneStr = ''.join(Known_Interactor_GeneList)
                Known_Interactors_Genes.append(Known_Interactor_GeneStr)
            print(''.join(candidateGene), '\t', CandidateENSG_data[1], '\t', CandidateENSG_data[2], '\t', CandidateENSG_data[3], '\t', len(CandidateENSG_data[5]), '\t', ','.join(Known_Interactors_Genes))
        else:
            # No known interactors indicated by '-' in the 3rd and fourth columns
            print(''.join(candidateGene), '\t', CandidateENSG_data[1], '\t', CandidateENSG_data[2], '\t', CandidateENSG_data[3], '\t', '-', '\t', '-')

    return

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
-----------------------------------------------------------------------------------------------------------------------
Program: Parses the patient Candidate Gene file(s) and Canonical Transcripts file. Using the Interactome file generated
         by Interactome.py, checks the number of interactors for each candidate gene. Next, checks the number of
         interactors that are known candidate genes and prints to STDOUT in .tsv format
-----------------------------------------------------------------------------------------------------------------------
The output consists of 6 columns in .tsv format:
 -> Candidate Gene
 -> Pathology Identifier
 -> Confidence score
 -> Number of Interactors with the candidate gene
 -> Number of Interactors that are known candidate genes
 -> A comma-seperated list of these known Interactors
-----------------------------------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inCandidateFile', metavar = "Input File", dest = "inCandidateFile", nargs = '+', help = 'Input File Name (Patient meta data file in .xlsx format)', required = True)
    required.add_argument('--inCanonicalFile', metavar = "Input File", dest = "inCanonicalFile", help = 'Canonical Transcripts file', required = True)
    required.add_argument('--inInteractome', metavar = "Input File", dest = "inInteractome", help = 'Input File Name (High-quality Interactome (.tsv) produced by Interactome.py)', required = True)

    args = file_parser.parse_args()
    Lead1_CandidateGene(args)

if __name__ == "__main__":
    # Logging to the file
    date = time.strftime("%Y_%m_%d-%H%M%S")
    Log_Format = "%(levelname)s %(asctime)s - %(message)s \n"
    logging.basicConfig(filename ='Lead-1_candidateGenes_%s.log' % date, filemode = 'a', format  = Log_Format, level = logging.DEBUG)
    main()
