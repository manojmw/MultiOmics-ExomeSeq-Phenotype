#!/usr/bin/python

import argparse, sys
import pandas as pd

###Function for creating transcripts_Gene dictionary###
# Takes tab-seperated canonical transcripts file as INPUT
# Creates a dictionary using 2 columns: ENSG and Gene
# Key -> Gene; Value - ENSG
# Returns the dictionary
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
        sys.exit("Missing required column title: 'ENSG' \n")
    elif not Gene_col >= 0:
        sys.exit("Missing required column title: 'GENE' \n")
    # else grabbed the required column indexes -> PROCEED

    # Parsing the Uniprot Primary Accession file
    for line in Canonical_File:
        line = line.rstrip('\n')
        CanonicalTranscripts_fields = line.split('\t')

        # Key -> ENSG
        # Value -> Gene
        ENSG_Gene_dict[CanonicalTranscripts_fields[Gene_col]] = CanonicalTranscripts_fields[ENSG_col]

    return ENSG_Gene_dict

###Function for extracting data from candidate gene file###
# Takes candidate gene file in .xlsx format as INPUT
# Extracts data from 3 columns - Gene, pathologyID and Confidence score
# Returns a DataFrame with the data from above 3 columns - Candidate_data
def CandidateGeneParser(inCandidateFile):

    # Input files
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

###Function for mapping Candidate Genes to canonical ENSGs###
# Calls the functions ENSG_Gene and CandidateGeneParser
# Parses the input files as described in the respective functions
# Maps the Genes to ENSGs
# Prints to STDOUT in .tsv format
# Output consists of 3 columns: ENSG, PathologyID and Confidence Score
def CandidateGene2ENSG(args):

    # Calling the function ENSG_Gene
    Transcripts_Gene_dict = ENSG_Gene(args.inCanonicalFile)

    # Calling the function CandidateGeneParser
    CandidateGene_data = CandidateGeneParser(args.inCandidateFile)

    lost_CandidateGene = 0

    for data in CandidateGene_data:
        if data[0] in Transcripts_Gene_dict.keys():
            print(Transcripts_Gene_dict.get(data[0]), '\t', data[1], '\t', data[2])
        else:
            lost_CandidateGene += 1

    print(lost_CandidateGene, file = sys.stderr)
    return


####Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
--------------------------------------------------------------------------------------------------------------------
Program: Parses the patient Candidate Gene file(s) and Canonical Transcripts file, processes it and prints to STDOUT
--------------------------------------------------------------------------------------------------------------------
The output consists of 3 columns in .tsv format:
 -> Gene name
 -> Pathology Identifier
 -> Confidence score
-----------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inCandidateFile', metavar = "Input File", dest = "inCandidateFile", nargs = '+', help = 'Input File Name (Patient meta data file in xlsx format)', required = True)
    required.add_argument('--inCanonicalFile', metavar = "Input File", dest = "inCanonicalFile", help = 'Canonical Transcripts file', required = True)

    args = file_parser.parse_args()
    CandidateGene2ENSG(args)

if __name__ == "__main__":
    main()
