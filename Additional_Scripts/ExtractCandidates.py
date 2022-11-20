#!/usr/bin/python

# manojmw
# 15 July 2022

import argparse
import logging
import sys

###########################################################

# Parses the patient samples file in .xlsx format
# Required columns are: 'Causal gene' & 'pathologyID' 
# (can be in any order, but they MUST exist)
#
# Returns a dictionary
# - Key: pathologyID
# - Value: list of causal genes assoicated with each pathology
def ExtractCandidates(args):

    logging.info("Starting to run...")
    
    # Input - list of filtered files containing potential candidates
    filtcandidate_files = args.infilteredcandidates

    # Intializing a list to store data associated with
    # each filtered candidate file 
    # The list consist of sublists
    # - Each sublist represents one file
    # - The first item of each sublist will be the name 
    # of the method associated with the file
    Candidates_Method_Data = [[] for i in range(len(filtcandidate_files))]

    # Initializing dictionary to store known diseases genes/candidates
    # key: Gene 
    # value: 1
    KnownCandidates_dict = {}

    # Data lines
    for i in range(len(filtcandidate_files)):
        
        # Get the pathology and method associated with the file
        # If the file doesn't exists in the present working directory
        # and a path to the file is provided instead (Ex: ~/Desktop/sammy/MMAF_1-hop_*.xlsx)
        if '/' in filtcandidate_files[i]:
            FName = filtcandidate_files[i].split('/')[-1]
            # get the pathology name associated with the current file
            current_patho = FName.split('_')[0]
            # get the identification method associated with the current file
            IMethod = FName.split('_')[1]
        else: # if the file is present in the present working directory 
            # get the pathology name associated with the current file
            current_patho = filtcandidate_files[i].split('_')[0]
            # get the identification method associated with the current file
            IMethod = filtcandidate_files[i].split('_')[1]
        
        # Store this data in the Candidates_Method_dict
        # value is a list which will later be filled with candidates
        Candidates_Method_Data[i].append(IMethod)
        
        filtcandidate_file = open(filtcandidate_files[i])

        filtcandidate_file_hLine = filtcandidate_file.readline() # Grabbing the header line

        filtcandidate_hFields = filtcandidate_file_hLine.split('\t')

        # Check the column headers and grab indexes of our columns of interest
        (Gene_index, KnownCandGene_index) = (-1,-1)

        for hI in range(len(filtcandidate_hFields)):
            if filtcandidate_hFields[hI] == 'GENE':
                Gene_index = hI
            elif filtcandidate_hFields[hI] == 'KNOWN_CANDIDATE_GENE':
                KnownCandGene_index = hI

        # sanity check
        if not Gene_index >= 0:
            logging.error("Missing required column title 'GENE' in the file: %s \n" % filtcandidate_files[i])
            sys.exit()
        elif not KnownCandGene_index >= 0:
            logging.error("Missing required column title 'KNOWN_CANDIDATE_GENE' in the file: %s \n" % filtcandidate_files[i])
            sys.exit()
        # else grabbed the required column indexes -> PROCEED
        
        # Data lines
        for line in filtcandidate_file:
            line = line.rstrip('\n')
            line_fields = line.split('\t')

            # store the candidates 
            Candidates_Method_Data[i].append(line_fields[Gene_index])

            # remove extra quotes
            line_fields[KnownCandGene_index] = line_fields[KnownCandGene_index].strip('"')

            # A gene can be associated with multiple pathologies
            # this will be present in the 'KNOWN_CANDIDATE_GENE' if the gene
            # is a known candidate
            # each patho will be seperated by a comma
            try: 
                Gene_pathologies = line_fields[KnownCandGene_index].split(',')
            except:
                Gene_pathologies = line_fields[KnownCandGene_index]
                
            # if a gene is known candidate, store it in KnownCandidates_dict dictionary
            if current_patho in Gene_pathologies:
                if not line_fields[Gene_index] in KnownCandidates_dict:
                    KnownCandidates_dict[line_fields[Gene_index]] = 1
        
        # Closing the file after parsing all the lines
        filtcandidate_file.close()
    
    # Initializing a dictionary to store Ver high and high confidence candidates
    HConfidence_Candidates_dict = {}
    
    # Classifiying candidates
    VHighConfidence_Candidates = sorted(list(set(Candidates_Method_Data[0]) & set(Candidates_Method_Data[1]) & set(Candidates_Method_Data[2])))
    print("# Very High Confidence_Candidates - %s" % [item[0] for item in Candidates_Method_Data])
    for gene in VHighConfidence_Candidates:
        HConfidence_Candidates_dict[gene] = 1
        if gene in KnownCandidates_dict:
            print(gene, "(Known Candidate Gene)")
        else:
            print(gene)

    for m1 in range(len(Candidates_Method_Data) - 1):
        for m2 in range(m1+1, len(Candidates_Method_Data)): 
            HighConfidence_Candidates = sorted(list(set(Candidates_Method_Data[m1]) & set(Candidates_Method_Data[m2])))
            print("\n# High Confidence_Candidates - %s" % [Candidates_Method_Data[m1][0], Candidates_Method_Data[m2][0]])
            for gene in HighConfidence_Candidates:
                if not gene in VHighConfidence_Candidates:
                    HConfidence_Candidates_dict[gene] = 1
                    if gene in KnownCandidates_dict:
                        print(gene, "(Known Candidate Gene)")
                    else:
                        print(gene)
    
    for CMData in Candidates_Method_Data:
        print("\n# Low Confidence_Candidates - [%s]" % CMData[0])
        for gene in CMData[1:]:
            if not gene in HConfidence_Candidates_dict:
                if gene in KnownCandidates_dict:
                    print(gene, "(Known Candidate Gene)")
                else:
                    print(gene)

    logging.info("All done, completed successfully!")

    return

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
-----------------------------------------------------------------------------------------------------------------------------
Program: Parses the filtered files containing candidates identified using three different methods (Ex: 1-hop approach, 
         clustering, etc). Prints to STDOUT the same list of candidates classified according to their confidence levels.
-----------------------------------------------------------------------------------------------------------------------------
-> The criteria for confidence level is described below:
    - Very high confidence candidates: Identified using all three methods 
    - High confidence candidates: Identified using any two methods 
    - Low confidence candidates: Identified using only a single method 

!!! Note !!!: 
- The file name should be in the format patho_method_*.tsv (Ex: MMAF_1-hop_220511.tsv, MMAF_R1DREAM_220511.tsv). This 
  is very important to identify the pathology & method the file is associated with. The name of the method should also
  be followed by an underscore ('_')

- The first 2 columns ('GENE' and 'KNOWN_CANDIDATE_GENE') must exist and is used to get the intersection and also to check
  if any of the candidates are already known disease genes.

- This script is written to identify the intersection of candidates identified using 3 different methods So, it is important 
  to provide all 3 result files for the script to work.

- All input files should be associated with the same pathology, else the script will produce wrong results while marking
  positive candidates.
-----------------------------------------------------------------------------------------------------------------------------

Arguments [defaults] -> Can be abbreviated to shortest unambiguous prefixes
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--infilteredcandidates', metavar = "Input File", dest = "infilteredcandidates", nargs = '+', help = 'Filtered files (.tsv) containing candidates', required=True)
 
    args = file_parser.parse_args()
    ExtractCandidates(args)

if __name__ == "__main__":
    # Logging to Standard Error
    Log_Format = "%(levelname)s - %(asctime)s - %(message)s \n"
    logging.basicConfig(stream = sys.stderr, format  = Log_Format, level = logging.DEBUG)
    main()
