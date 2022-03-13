#!/usr/bin/python

import argparse, sys
import pandas as pd
import scipy.stats as stats
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
#
# Maps Gene name to ENSG
# Returns 2 lists:
# The first list has 3 items (in each sublist:
# - ENSG
# - pathologyID
# - Confidence score)
# The second list contains all the pathologies
def CandidateGene2ENSG(ENSG_Gene_dict, CandidateGene_data):

    # List of pathologies
    pathologies_list = []

    # Initializing output list
    candidateENSG_out_list = []

    # Counter for canidate genes not found in the canonical transcripts file
    lost_CandidateGene = 0

    for data in CandidateGene_data:
        # data[0] -> Candidate Gene
        # data[1] -> pathologyID
        # data[2] -> Confidence Score
        if data[0] in ENSG_Gene_dict.keys(): # Get the corresponding ENSG
            candidateENSG_out = [ENSG_Gene_dict.get(data[0]),  data[1],  data[2]]
            candidateENSG_out_list.append(candidateENSG_out)
        else:
            lost_CandidateGene += 1

        if not data[1] in pathologies_list:
            pathologies_list.append(data[1])

    logging.debug("\nNo. of candidate genes not found in the Canonical transcripts file: %d" % lost_CandidateGene)

    return candidateENSG_out_list, pathologies_list

###########################################################

# Parses the candidateENSG_out_list & pathologies_list
# Counts the total number of candidate genes
# associated with each pathology
# Returns a list with the total count (for each pathology)
def CountCandidateGenes(candidateENSG_out_list, pathologies_list):

    # List for counting total candidate genes
    # associated with each pathology
    pathology_CandidateCount = [0] * len(pathologies_list)

    for candidateGenedata in candidateENSG_out_list:
        for i in range(len(pathologies_list)):
            if candidateGenedata[1] == pathologies_list[i]:
                pathology_CandidateCount[i] += 1

    return pathology_CandidateCount

###########################################################

# Parses the High-quality Interactome produced by Interactome.py
# Required columns:
# First column - ENSG of Protein A
# Second column - ENSG of Protein B
#
# Returns 2 lists:
# First list contains 2 items (in each sub list):
# - ENSG of Protein A
# - ENSG of Protein B
# Second list contains all the interacting proteins from the interactome
def Interacting_Proteins(inInteractome):

    # Input - Interactome file
    Interactome_File = open(inInteractome)

    # Initializing Interactome dictionary
    Interactome_list = []

    # List of all interactings from the Interactome
    All_Interactors_list = []

    for line in Interactome_File:
        line = line.rstrip('\n')

        Interactome_fields = line.split('\t')

        Interacting_Proteins = [Interactome_fields[0], Interactome_fields[1]]

        Interactome_list.append(Interacting_Proteins)

        # Storing all the interactors in All_Interactors_list
        if not Interactome_fields[0] in All_Interactors_list:
            All_Interactors_list.append(Interactome_fields[0])
        elif not Interactome_fields[1] in All_Interactors_list:
            All_Interactors_list.append(Interactome_fields[1])

    # Closing the file
    Interactome_File.close()

    return Interactome_list, All_Interactors_list

###########################################################

# Parses the Interactome_list & All_Interactors_list returned
# by the function: Interacting_Proteins
# Checks the number of interactors for each gene
# Checks the number of known interactors
# using the candidateGene_out_list returned by the function: CandidateGene2ENSG
#
# Prints to STDOUT in .tsv format
# The output consists of sublist(s) per each line
# Each sublist represents one pathology/phenotype and consists of:
# - Gene
# - Number of Interactors
# - No. of Known interactors
# - Benjamini-Hochberg adjusted P-value
def Interactors_PValue(args):

    # Calling the functions
    CandidateGene_data = CandidateGeneParser(args.inCandidateFile)
    ENSG_Gene_dict = ENSG_Gene(args.inCanonicalFile)
    (Interactome_list, All_Interactors_list) = Interacting_Proteins(args.inInteractome)
    (candidateENSG_out_list, pathologies_list) = CandidateGene2ENSG(ENSG_Gene_dict, CandidateGene_data)
    pathology_CandidateCount = CountCandidateGenes(candidateENSG_out_list, pathologies_list)

    total_human_ENSG = 22000 # Approximate count, will be corrected later using Uniprot file

    # Initializing first output list with p-values
    # each sublist represents one pathology
    patho_p_value = []

    for i in range(len(pathologies_list)):

        # Initializing a list to store data for each pathology
        Output_eachPatho = []

        # Checking the number of interactors for each gene
        for ENSG in All_Interactors_list:

            # List of interactors
            Interactors = []

            # List for known interactor
            Known_Interactors = 0

            for Proteins in Interactome_list:
                # If Protein_A is the first protein
                if (ENSG == Proteins[0]):
                    # Get the interacting protein
                    if not Proteins[1] in Interactors:
                        Interactors.append(Proteins[1])
                # If Protein_A is the Second protein
                elif (ENSG == Proteins[1]):
                    if not Proteins[0] in Interactors:
                        # Get the interacting protein
                        Interactors.append(Proteins[0])

            # Checking if the interactor is a known ENSG (candidate ENSG)
            for interactor in Interactors:
                for candidateENSG in candidateENSG_out_list:
                    if interactor in candidateENSG:
                        if candidateENSG[1] == pathologies_list[i]:
                            Known_Interactors += 1

            # Applying Fisher's exact test to calculate p-values
            raw_data = [[Known_Interactors, len(Interactors)],[pathology_CandidateCount[i], total_human_ENSG]]
            (odd_ratio, p_value) = stats.fisher_exact(raw_data)

            Output_eachENSG = [ENSG, len(Interactors), Known_Interactors, p_value]

            Output_eachPatho.append(Output_eachENSG)

        patho_p_value.append(Output_eachPatho)

    # Sorting the p-values for each pathology
    for i in range(len(patho_p_value)):
        patho_p_value[i].sort(key = lambda x:x[3])

        # Calculating Benjamini-Hochberg corrected p-values
        for data in patho_p_value[i]:

            # Rank of p-value -> patho_p_value[i].index(data) + 1
            # (+1 because list index starts from 0)

            # p-value -> data[3]

            # Total number of tests -> len(patho_p_value[i])

            # Benjamini Hochberg corrected p-value = p-value*Total number of tests/Rank of p-value

            BH_p_value = round((data[3] * len(patho_p_value[i]))/(patho_p_value[i].index(data)+1), 2)
            data[3] = BH_p_value

    # Printing the header for the output
    print('\t'.join(pathology for pathology in pathologies_list))

    for i in range(len(patho_p_value)):
        # * allows us to unpack the sublists
        # So that we can access indivdual elements of the sublist
        # zip to get a list of the first element from the sublists
        # list because zip in python 3 doesn't return a list
        # Once we get a list, we iterate and add a tab seperator
        out_tuple = (list(zip(*patho_p_value))[i])
        out_tabsep = '\t'.join(str(value) for value in out_tuple)
        print(out_tabsep)

    return

# ###########################################################
#
# # Parses the list (candGene_Interactors_list)
# # returned by the function: Lead1_CandidateENSG
# # Maps the ENSGs in the list to Gene names using the dictionary: ENSG_Gene_dict
# #
# # Prints to STDOUT in .tsv format
# # Output consists of 6 columns:
# # - Candidate Gene
# # - Pathology Identifier
# # - Confidence Score
# # - Count of Proteins interacting with Candidate Gene
# # - Count of Interacting genes that are known candidate Genes
# # - List of Known Interacting genes
# def Lead1_CandidateGene(args):
#
#     candENSG_Interactors_list = Lead1_CandidateENSG(args.inCanonicalFile, args.inCandidateFile, args.inInteractome)
#     ENSG_Gene_dict = ENSG_Gene(args.inCanonicalFile)
#
#     # Printing the header for the output
#     header = ['Candidate_Gene', 'pathologyID', 'Confidence_Score', 'No_of_Interactors', 'No_of_KnownInteractors', 'List_KnownInteractors']
#     print('\t'.join(header))
#
#     for CandidateENSG_data in candENSG_Interactors_list:
#
#         # known interactor genes
#         Known_Interactors_Genes = []
#
#         candidateGene = [GeneName for (GeneName, ENSG) in ENSG_Gene_dict.items() if ENSG == CandidateENSG_data[0]]
#
#         # If there are known interactors in the candGene_Interactors_list
#         # Get the Gene name of the known interactor using the dictionary (ENSG_Gene_dict)
#         for Interactors_ENSG in CandidateENSG_data[5]:
#             Known_Interactor_GeneList = [GeneName for (GeneName, ENSG) in ENSG_Gene_dict.items() if ENSG == Interactors_ENSG]
#             Known_Interactor_GeneStr = ''.join(Known_Interactor_GeneList)
#             Known_Interactors_Genes.append(Known_Interactor_GeneStr)
#
#         # Print to STDOUT
#         print(''.join(candidateGene), '\t', CandidateENSG_data[1], '\t', CandidateENSG_data[2], '\t', CandidateENSG_data[3], '\t', len(CandidateENSG_data[5]), '\t', ','.join(Known_Interactors_Genes))
#
#     return

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
-----------------------------------------------------------------------------------------------------------------------
Program: Parses the Interactome file generated by Interactome.py, checks the number of interactors for each gene. Next,
         parses the patient Candidate Gene file(s) and Canonical Transcripts file. Checks the number of interactors
         that are known candidate genes, calculates P-values and prints to STDOUT in .tsv format
-----------------------------------------------------------------------------------------------------------------------
The output consists of sublist(s) per each line. Each sublist represents one pathology/phenotype and consists of:
 -> Gene
 -> Number of Interactors
 -> No. of Known interactors
 -> Benjamini-Hochberg adjusted P-value
-----------------------------------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inCandidateFile', metavar = "Input File", dest = "inCandidateFile", nargs = '+', help = 'Input File Name (Patient meta data file in .xlsx format)', required = True)
    required.add_argument('--inCanonicalFile', metavar = "Input File", dest = "inCanonicalFile", help = 'Canonical Transcripts file', required = True)
    required.add_argument('--inInteractome', metavar = "Input File", dest = "inInteractome", help = 'Input File Name (High-quality Interactome (.tsv) produced by Interactome.py)', required = True)

    args = file_parser.parse_args()
    Interactors_PValue(args)

if __name__ == "__main__":
    # Logging to the file
    date = time.strftime("%Y_%m_%d-%H%M%S")
    Log_Format = "%(levelname)s %(asctime)s - %(message)s \n"
    logging.basicConfig(filename ='Lead-1_candidateGenes_%s.log' % date, filemode = 'a', format  = Log_Format, level = logging.DEBUG)
    main()
