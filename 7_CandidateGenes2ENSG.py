#!/usr/bin/python

import argparse, sys
import pandas as pd
import logging
import time

###########################################################

# Parses tab-seperated canonical transcripts file
# Required columns are: 'ENSG' and 'GENE' (can be in any order,
# but they MUST exist)
# Returns a dictionary:
# Key -> ENSG; Value -> Gene
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

        # Key -> ENSG
        # Value -> Gene
        ENSG_Gene_dict[CanonicalTranscripts_fields[Gene_col]] = CanonicalTranscripts_fields[ENSG_col]

    return ENSG_Gene_dict

###########################################################

###Function for extracting data from candidate gene file###
# Takes candidate gene file in .xlsx format as INPUT
# Extracts data from 3 columns - Gene, pathologyID and Confidence score
# Returns a DataFrame with the data from above 3 columns - Candidate_data
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

###Function for mapping Candidate Genes to canonical ENSGs###
# Calls the functions ENSG_Gene and CandidateGeneParser
# Parses the input files as described in the respective functions
# Maps the Genes to ENSGs
# Prints to STDOUT in .tsv format
# Output consists of 3 columns: ENSG, PathologyID and Confidence Score
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

    return Interactome_list

###########################################################

# Parses the canidateGene_out_list returned by the function: CandidateGene2ENSG
# Checks the Number of interactors for each candidate gene
# using the Interactome_list returned by the function: Interacting_Proteins
# Checks the number of Interactors that are known candidate genes
# Prints to STDOUT in .tsv format
# Output consists of 5 columns:
# - Candidate Gene
# - Count of Proteins interacting with Candidate Gene
# - List of ENSGs of Interacting proteins
# - Count of Interacting proteins that are known candidate Genes
# - List of Candidate Interacting Proteins
def Lead1_CandidateGenes(args):

    # Calling the functions
    canidateGene_out_list = CandidateGene2ENSG(args.inCanonicalFile, args.inCandidateFile)
    Interactome_list = Interacting_Proteins(args.inInteractome)

    # Dictionary for storing the Candidate Genes and Interactors
    candGene_Interactors_list = []


    for candidateGene in canidateGene_out_list:

        Interactors = []

        for Proteins in Interactome_list:
            if (candidateGene[0] == Proteins[0]):
                Interactors.append(Proteins[1])
                candGene_Interactors = [candidateGene[0], candidateGene[1], len(Interactors), ''.join(Interactors)]
                candGene_Interactors_list.append(candGene_Interactors)
            elif (candidateGene[0] == Proteins[1]):
                Interactors.append(Proteins[0])
                candGene_Interactors = [candidateGene[0], candidateGene[1], len(Interactors), ','.join(Interactors)]
                candGene_Interactors_list.append(candGene_Interactors)


    return print(candGene_Interactors_list)



###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
--------------------------------------------------------------------------------------------------------------------
Program: Parses the patient Candidate Gene file(s) and Canonical Transcripts file, processes it and prints to STDOUT
--------------------------------------------------------------------------------------------------------------------
The output consists of 3 columns in .tsv format:
 -> Candidate ENSG
 -> Pathology Identifier
 -> Confidence score
--------------------------------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inCandidateFile', metavar = "Input File", dest = "inCandidateFile", nargs = '+', help = 'Input File Name (Patient meta data file in xlsx format)', required = True)
    required.add_argument('--inCanonicalFile', metavar = "Input File", dest = "inCanonicalFile", help = 'Canonical Transcripts file', required = True)
    required.add_argument('--inInteractome', metavar = "Input File", dest = "inInteractome", help = 'Input File Name (High-quality Interactome (.tsv) produced by Interactome.py)', required = True)

    args = file_parser.parse_args()
    Lead1_CandidateGenes(args)

if __name__ == "__main__":
    # Logging to the file
    date = time.strftime("%Y_%m_%d-%H%M%S")
    Log_Format = "%(levelname)s %(asctime)s - %(message)s \n"
    logging.basicConfig(filename ='Lead-1_candidateGenes_%s.log' % date, filemode = 'a', format  = Log_Format, level = logging.DEBUG)
    main()
