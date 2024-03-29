#!/usr/local/bin/python3

# manojmw
# 18 Nov 2022

import sys, argparse
import logging
import gzip
import re

###########################################################

# Parses tab-seperated canonical transcripts file (.gz or non .gz)
# Required columns are: 'ENSG' and 'GENE' (can be in any order,
# but they MUST exist)
#
# Returns a dictionary:
# - Key -> Gene
# - Value -> ENSG
def Gene_ENSG(incanonical):

    logging.info("Starting to run...")

    # Dictionary for storing ENSG
    # and Gene data
    Gene_ENSG_dict = {}

    # Opening canonical transcript file (gzip or non-gzip)
    try:
        if incanonical.endswith('.gz'):
            Canonical_File = gzip.open(incanonical, 'rt')
        else:
            Canonical_File = open(incanonical)
    except IOError:
        sys.exit("Error: Failed to read the Canonical transcript file: %s" % incanonical)

    # Grabbing the header line
    Canonical_header_line = Canonical_File.readline()

    Canonical_header_fields = Canonical_header_line.split('\t')

    # Check the column headers and grab indexes of our columns of interest
    (ENSG_index, Gene_index) = (-1,-1)

    for i in range(len(Canonical_header_fields)):
        if Canonical_header_fields[i] == 'ENSG':
            ENSG_index = i
        elif Canonical_header_fields[i] == 'GENE':
            Gene_index = i

    if not ENSG_index >= 0:
        sys.exit("Error: Missing required column title 'ENSG' in the file: %s \n" % incanonical)
    elif not Gene_index >= 0:
        sys.exit("Error: Missing required column title 'GENE' in the file: %s \n" % incanonical)
    # else grabbed the required column indexes -> PROCEED

    # Data lines
    for line in Canonical_File:
        line = line.rstrip('\n')
        CanonicalTranscripts_fields = line.split('\t')

        # Key -> ENSG
        # Value -> Gene

        Gene_ENSG_dict[CanonicalTranscripts_fields[Gene_index]] = CanonicalTranscripts_fields[ENSG_index]

    # Closing the file
    Canonical_File.close()

    return Gene_ENSG_dict

###########################################################

# Parses the tab-seperated UniProt Primary Accession file
# produced by 1_Uniprot_parser.py
#
# Required columns are: 'Alt_GeneID','ENSGs', 'GeneName'
# and 'Function' (can be in any order,
# but they MUST exist)

# Returns two dictionaries:
# First dictionary - ENSG_AltGID_dict
# - Key: ENSG
# - Value: Alternate Gene ID
# Second dictionary - AltGID_GeneNFunc_dict
# - Key: Alternate Gene ID
# - Value: List(s) containing Gene name and Function

def Uniprot_parse(inuniprot):

    # Dictionary for storing ENSG (key)
    # and Alternate Gene ID (value) (Ex: HGNC ID, MGI ID, etc)
    ENSG_AltGID_dict = {}

    # Dictionary for storing Alternate Gene ID (key)
    # and value (Gene name and Function)
    AltGID_GeneNFunc_dict = {}

    Uniprot_File = open(inuniprot)

    # Grabbing the header line
    Uniprot_header = Uniprot_File.readline()

    Uniprot_header = Uniprot_header.rstrip('\n')

    Uniprot_header_fields = Uniprot_header.split('\t')

    # Check the column header and grab indexes of our columns of interest
    (AltGeneID_index, ENSG_index, GeneName_index, Function_index) = (-1, -1, -1, -1)

    for i in range(len(Uniprot_header_fields)):
        if Uniprot_header_fields[i] == 'Alt_GeneID':
            AltGeneID_index = i
        elif Uniprot_header_fields[i] == 'ENSGs':
            ENSG_index = i
        elif Uniprot_header_fields[i] == 'GeneName':
            GeneName_index = i
        elif Uniprot_header_fields[i] == 'Function':
            Function_index = i

    # Sanity check
    if not AltGeneID_index >= 0:
        sys.exit("Error: Missing required column title 'Alt_GeneID' in the file: %s \n" % inuniprot)
    elif not ENSG_index >= 0:
        sys.exit("Error: Missing required column title 'ENSGs' in the file: %s \n" % inuniprot)
    elif not GeneName_index >= 0:
        sys.exit("Error: Missing required column title 'GeneName' in the file: %s \n" % inuniprot)
    elif not Function_index >= 0:
        sys.exit("Error: Missing required column title 'Function' in the file: %s \n" % inuniprot)
    # else grabbed the required column indices -> PROCEED

    # Data lines
    for line in Uniprot_File:
        line = line.rstrip('\n')
        Uniprot_fields = line.split('\t')

        # ENSG column  - This is a single string containing comma-seperated ENSGs
        # So we split it into a list that can be accessed later
        UniProt_ENSGs = Uniprot_fields[ENSG_index].split(',')
        
        if len(UniProt_ENSGs) == 1:
            # ENSG_HGNC_dict:
            # - key: ENSG
            # - value: AltGeneID
            ENSG_AltGID_dict[Uniprot_fields[ENSG_index]] = Uniprot_fields[AltGeneID_index]
        elif len(UniProt_ENSGs) > 1:
            for ENSG in UniProt_ENSGs:
                ENSG_AltGID_dict[ENSG] = Uniprot_fields[AltGeneID_index]
        
        Alt_GeneID = Uniprot_fields[AltGeneID_index]
        GeneName = Uniprot_fields[GeneName_index]
        Function = Uniprot_fields[Function_index]
        # Storing Alt_GeneID, Gene name and function data in AltGID_GeneNFunc_dict
        if Alt_GeneID in AltGID_GeneNFunc_dict:
            # Alt_GeneID already exists, append the new sublist to the values
            AltGID_GeneNFunc_dict[Alt_GeneID].append([GeneName, Function])
        else: 
             AltGID_GeneNFunc_dict[Alt_GeneID] = [[GeneName, Function]]
        
    # Closing the file
    Uniprot_File.close()

    return ENSG_AltGID_dict, AltGID_GeneNFunc_dict

###########################################################

# Parses Candidate Results file produced by ExtractCandidates.py
# Also parses dictionaries (Gene_ENSG_dict, ENSG_AltGID_dict)
# returned by Uniprot_parse function
#
# Returns a dictionary:
# - Key: Name of the Confidence level
# - Value: A containing:
#        - Key: Gene
#        - Value

def Parse_CandidateResults(incandidateresult, Gene_ENSG_dict, ENSG_AltGID_dict):

    CandResult_File = open(incandidateresult)

    # Dictionary to store Candidate results
    # - key: Confidence level name
    # - value: Dictionary containing:
    #          - Key: Gene; 
    #          - value: 1 or "(Known Candidate Gene)" if
    #                   the gene is a known candidate
    CandResult_dict = {}

    # In candidate results file, the data is of the format:
    # # Very High Confidence_Candidates
    # Gene 1
    # Gene 2
    #
    # # High Confidence_Candidates
    # Gene 3
    # Gene 4
    #
    # An emty line indicates the end of candidate genes in
    # the given confidence level
    # inConfidence == True; when we are in a Confidence level block, 
    # False otherwise; inConfidence == False (when we see an empty line
    # - indicates end of candidate genes in the given confidence level)
    inConfidence = False

    # Data lines
    for line in CandResult_File:
        if line.startswith("#"):
            inConfidence = True
            ConfidenceL_key = line.rstrip('\n')
            CandResult_dict[ConfidenceL_key] = {}
        elif inConfidence:
            line = line.rstrip('\r\n')
            if re.match(r'(\S+)', line):
                Candidate_Gene = re.match(r'(\S+)', line).group(1)
                # print(Candidate_Gene)
                if Candidate_Gene in Gene_ENSG_dict:
                    # Get ENSG ID for the candidate Gene
                    ENSG = Gene_ENSG_dict[Candidate_Gene]
                    # Get HGNC ID for the ENSG
                    if ENSG in ENSG_AltGID_dict:
                        HGNC_ID = ENSG_AltGID_dict[ENSG]
                        # If the gene is a known candidate, then store the value
                        # as "(Known Candidate Gene)"
                        if "(Known Candidate Gene)" in line:
                             CandResult_dict[ConfidenceL_key][HGNC_ID] = "(Known Candidate Gene)"
                        else: # If the gene is not a known candidate
                            CandResult_dict[ConfidenceL_key][HGNC_ID] = 1
                    else:
                        logging.info("Could not get HGNC ID for the candidate Gene: " + line + 'with ENSG ID: ' + ENSG)
                else:
                    logging.info("Could not get ENSG ID for the candidate Gene: " + line)
                
        elif (line == " "): # we have reached the end of the current confidence level
            inConfidence = False
            continue
    
    return CandResult_dict

###########################################################

# Parse the tab-separated Orthology file (.gz or non .gz)
# retrieved from Alliance of Genome Resources
#
# Returns 2 dictionaries:
# First Dictionary (Gene1_dict) contains:
# - Key: Gene1ID
# - Value: List contining following data for both Gene 1 and Gene 2:
    #          - Symbol
    #          - TaxID
    #          - Algorithm Match count
    #          - Total number of algorithms
    #          - Best Score
    #          - Best Rev Score
    #          - Additionally contains Gene2ID
# Second dictionary (Gene2_dict) contains:
# - Key: Gene2ID
# - Value: same data as Gene1_dict, but additionally
#          contains Gene1ID istead of Gene2ID
def Parse_OrthologFile(inortholog):

    # Opening orthology file (gzip or non-gzip)
    try:
        if inortholog.endswith('.gz'):
            OrthoF = gzip.open(inortholog, 'rt')
        else:
            OrthoF = open(inortholog)
    except IOError:
        sys.exit("Error: Failed to read the Orthology file: %s" % inortholog)

    # Dictionaries to store orthology data
    Gene1_dict = {}
    Gene2_dict = {}

    # Column indices
    (Gene1ID_index, Gene1Symbol_index, Gene1TaxID_index, Gene2ID_index, Gene2Symbol_index, Gene2TaxID_index,
    AlgoMatchC_index, TotalAlgoC_index, BestScore_index, BestRevScore_index) = (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1)

    for line in OrthoF:

        line = line.rstrip('\r\n')

        if line.startswith('#'):
            continue # Skip metadata lines

        elif line.startswith('Gene1ID'): # Header line
            Ortho_header_fields = line.split('\t')

            for i in range(len(Ortho_header_fields)):
                if Ortho_header_fields[i] == 'Gene1ID':
                    Gene1ID_index = i
                elif Ortho_header_fields[i] == 'Gene1Symbol':
                    Gene1Symbol_index = i
                elif Ortho_header_fields[i] == 'Gene1SpeciesTaxonID':
                    Gene1TaxID_index = i
                elif Ortho_header_fields[i] == 'Gene2ID':
                    Gene2ID_index = i
                elif Ortho_header_fields[i] == 'Gene2Symbol':
                    Gene2Symbol_index = i
                elif Ortho_header_fields[i] == 'Gene2SpeciesTaxonID':
                    Gene2TaxID_index = i
                elif Ortho_header_fields[i] == 'AlgorithmsMatch':
                    AlgoMatchC_index = i
                elif Ortho_header_fields[i] == 'OutOfAlgorithms':
                    TotalAlgoC_index = i
                elif Ortho_header_fields[i] == 'IsBestScore':
                    BestScore_index = i
                elif Ortho_header_fields[i] == 'IsBestRevScore':
                    BestRevScore_index = i

            # Sanity check
            if not Gene1ID_index >= 0:
                sys.exit("Error: Missing required column title 'Gene1ID' in the file: %s \n" % inortholog)
            elif not Gene1Symbol_index >= 0:
                sys.exit("Error: Missing required column title 'Gene1Symbol' in the file: %s \n" % inortholog)
            elif not Gene1TaxID_index >= 0:
                sys.exit("Error: Missing required column title 'Gene1SpeciesTaxonID' in the file: %s \n" % inortholog)
            elif not Gene2ID_index >= 0:
                sys.exit("Error: Missing required column title 'Gene2ID' in the file: %s \n" % inortholog)
            elif not Gene2Symbol_index >= 0:
                sys.exit("Error: Missing required column title 'Gene2Symbol' in the file: %s \n" % inortholog)
            elif not Gene2TaxID_index >= 0:
                sys.exit("Error: Missing required column title 'Gene2SpeciesTaxonID' in the file: %s \n" % inortholog)
            elif not AlgoMatchC_index >= 0:
                sys.exit("Error: Missing required column title 'AlgorithmsMatch' in the file: %s \n" % inortholog)
            elif not TotalAlgoC_index >= 0:
                sys.exit("Error: Missing required column title 'OutOfAlgorithms' in the file: %s \n" % inortholog)
            elif not BestScore_index >= 0:
                sys.exit("Error: Missing required column title 'IsBestScore' in the file: %s \n" % inortholog)
            elif not BestRevScore_index >= 0:
                sys.exit("Error: Missing required column title 'IsBestRevScore' in the file: %s \n" % inortholog)
            # else grabbed the required column indices -> PROCEED
        
        # Data Lines
        else:
            orthodata_fields = line.split('\t')

            Gene1ID = orthodata_fields[Gene1ID_index]
            Gene1Symbol = orthodata_fields[Gene1Symbol_index]
            Gene1SpeciesTaxonID = orthodata_fields[Gene1TaxID_index]
            Gene2ID = orthodata_fields[Gene2ID_index]
            Gene2Symbol = orthodata_fields[Gene2Symbol_index]
            Gene2SpeciesTaxonID = orthodata_fields[Gene2TaxID_index]
            AlgorithmsMatch = orthodata_fields[AlgoMatchC_index]
            OutOfAlgorithms = orthodata_fields[TotalAlgoC_index]
            IsBestScore = orthodata_fields[BestScore_index]
            IsBestRevScore = orthodata_fields[BestRevScore_index]

            # Here, the list of candidates in the candidate result file are human genes (NCBITaxon:9606) 
            # and we are interested in finding the ortholog for Yeast Gene (NCByoITaxon:559292)
            # So, you will need to modify the script based on the model organism you are interested in
            # NOTE: you will also need to modify the 1_Uniprot_parser.py script as well for this to work
            # You can find the instructions for that in the 1_Uniprot_parser.py script

            if (Gene1SpeciesTaxonID == 'NCBITaxon:9606') and (Gene2SpeciesTaxonID == 'NCBITaxon:559292'):
                if not Gene1ID in Gene1_dict:
                    Gene1_dict[Gene1ID] = [[Gene1Symbol, Gene1SpeciesTaxonID, Gene2Symbol, Gene2SpeciesTaxonID, 
                                        (AlgorithmsMatch, OutOfAlgorithms, IsBestScore, IsBestRevScore), Gene2ID]]
                else:    
                    Gene1_dict[Gene1ID].append([Gene1Symbol, Gene1SpeciesTaxonID, Gene2Symbol, Gene2SpeciesTaxonID, 
                                            (AlgorithmsMatch, OutOfAlgorithms, IsBestScore, IsBestRevScore), Gene2ID])
            
            if (Gene1SpeciesTaxonID == 'NCBITaxon:559292') and (Gene2SpeciesTaxonID == 'NCBITaxon:9606'):
                if not Gene2ID in Gene2_dict:
                    Gene2_dict[Gene2ID] = [[Gene1Symbol, Gene1SpeciesTaxonID, Gene2Symbol, Gene2SpeciesTaxonID, 
                                        (AlgorithmsMatch, OutOfAlgorithms, IsBestScore, IsBestRevScore), Gene1ID]]
                else:    
                    Gene2_dict[Gene2ID].append([Gene1Symbol, Gene1SpeciesTaxonID, Gene2Symbol, Gene2SpeciesTaxonID, 
                                          (AlgorithmsMatch, OutOfAlgorithms, IsBestScore, IsBestRevScore), Gene1ID])       
                                
    return Gene1_dict, Gene2_dict

###########################################################

# Parses the dictionaries returned by above functions
# and searches for model organism orthologs for the Human candidate genes
# 
# Prints to STDOUT in .tsv format
# The line with name of the confidence level starts with an '#'
# The following lines include the Human candidate gene along
# with the following information:
# - Name of the Human candidate gene (including Synonyms)
# - HGNC ID
# - Whether the candidate gene is already a known candidate
#   for the phenotype interested in
# - Human candidate gene function
# - Model organism ortholog (if available)
# - Model organism Gene ID
# - Orthology match scores
# - Model organsim Gene function
def FindOrthologs(args):

    # Calling the functions
    Gene_ENSG_dict = Gene_ENSG(args.incanonical)
    (ENSG_AltGID_dict, AltGID_GeneNFunc_dict) = Uniprot_parse(args.inuniprot)
    CandResult_dict = Parse_CandidateResults(args.incandidateresult, Gene_ENSG_dict, ENSG_AltGID_dict)
    (Gene1_dict, Gene2_dict) = Parse_OrthologFile(args.inortholog)  

    header = ["HumanGene", "HGNC_ID", "Known_Candidate", "HumanGene_Function", "YeastOrtholog", "SGD_ID", "Orthology_MatchScores", "YeastGene_Function"]
    print('\t'.join(header))

    # NOTE: Again, the below code needs to modified based on the model organism you are interested in (particularly TaxIDs)

    # Finding the orthologs for our candidate genes
    for ConfidenceL_key in CandResult_dict:
        print(ConfidenceL_key)
        for Candidate_Gene in CandResult_dict[ConfidenceL_key]:
            if Candidate_Gene in Gene1_dict:
                for data in Gene1_dict[Candidate_Gene]:
                    
                    if data[1] != 'NCBITaxon:9606':
                        logging.error("Candidate gene is a human gene in both candidateresult file and orthology file, but the tax ID is not '9606', Impossible!!!")
                        sys.exit()

                    if data[3] == 'NCBITaxon:559292':
                        output = []
                        if Candidate_Gene in AltGID_GeneNFunc_dict:
                            for CandGene_data in AltGID_GeneNFunc_dict[Candidate_Gene]:
                                output.append(CandGene_data[0]) # append gene name from uniprot
                                output.append(Candidate_Gene) # append HGNC ID
                                # Check if the candidate gene is a known candidate
                                if CandResult_dict[ConfidenceL_key][Candidate_Gene] == "(Known Candidate Gene)":
                                    output.append("YES")
                                else:
                                    output.append(" ")   
                                output.append(CandGene_data[1]) # append function (if available)
                        else:
                            output.append(data[0]) # append gene name from orthology file
                            output.append(Candidate_Gene) # append HGNC ID
                            # Check if the candidate gene is a known candidate
                            if CandResult_dict[ConfidenceL_key][Candidate_Gene] == "(Known Candidate Gene)":
                                    output.append("YES")
                            else:
                                output.append(" ") 
                            output.append(" ") # no function
                        if data[-1] in AltGID_GeneNFunc_dict:
                            for OrthoGene_data in AltGID_GeneNFunc_dict[data[-1]]:
                                output.append(OrthoGene_data[0]) # append gene name from uniprot
                                output.append(data[-1] ) # append SGD ID
                                output.append(str(data[4])) # append ortholog match details
                                output.append(OrthoGene_data[1]) # append function (if available)
                        else:
                            output.append(data[2]) # append gene name from orthology file
                            output.append(data[-1] ) # append SGD ID
                            output.append(str(data[4])) # append ortholog match details
                            output.append(" ") # no function

                        print("\t".join(output))
                
            if Candidate_Gene in Gene2_dict:
                for data in Gene2_dict[Candidate_Gene]:

                    if data[3] != 'NCBITaxon:9606':
                        logging.error("Candidate gene is a human gene in both candidateresult file and orthology file, but the tax ID is not '9606', Impossible!!!")
                        sys.exit()

                    if data[1] == 'NCBITaxon:559292':
                        output = []
                        if Candidate_Gene in Gene1_dict:
                            for Gene1_dict_data in Gene1_dict[Candidate_Gene]:
                                if Gene1_dict_data[4] == data[4]: # same match scores data           
                                    continue  # avoid adding redundant data
                        else:
                            if Candidate_Gene in AltGID_GeneNFunc_dict:
                                for CandGene_data in AltGID_GeneNFunc_dict[Candidate_Gene]:
                                    output.append(CandGene_data[0]) # append gene name from uniprot
                                    output.append(Candidate_Gene) # append HGNC ID
                                    # Check if the candidate gene is a known candidate
                                    if CandResult_dict[ConfidenceL_key][Candidate_Gene] == "(Known Candidate Gene)":
                                        output.append("YES")
                                    else:
                                        output.append(" ") 
                                    output.append(CandGene_data[1]) # append function (if available)
                            else:
                                output.append(data[2]) # append gene name from orthology file
                                output.append(Candidate_Gene) # append HGNC ID
                                # Check if the candidate gene is a known candidate
                                if CandResult_dict[ConfidenceL_key][Candidate_Gene] == "(Known Candidate Gene)":
                                    output.append("YES")
                                else:
                                    output.append(" ") 
                                output.append(" ") # no function
                            if data[-1] in AltGID_GeneNFunc_dict:
                                for OrthoGene_data in AltGID_GeneNFunc_dict[data[-1]]:
                                    output.append(OrthoGene_data[0]) # append gene name from uniprot
                                    output.append(data[-1] ) # append SGD ID
                                    output.append(str(data[4])) # append ortholog match details
                                    output.append(OrthoGene_data[1]) # append function (if available)
                            else:
                                output.append(data[0]) # append gene name from orthology file
                                output.append(data[-1] ) # append SGD ID
                                output.append(str(data[4])) # append ortholog match details
                                output.append(" ") # no function

                            print("\t".join(output))

            if (Candidate_Gene not in Gene1_dict) and (Candidate_Gene not in Gene2_dict):
                output = []
                if Candidate_Gene in AltGID_GeneNFunc_dict:
                    for CandGene_data in AltGID_GeneNFunc_dict[Candidate_Gene]:
                        output.append(CandGene_data[0]) # append gene name from uniprot
                        output.append(Candidate_Gene) # append HGNC ID
                        # Check if the candidate gene is a known candidate
                        if CandResult_dict[ConfidenceL_key][Candidate_Gene] == "(Known Candidate Gene)":
                            output.append("YES")
                        else:
                            output.append(" ") 
                        output.append(CandGene_data[1]) # append function (if available)
                        output = output + [" ", " ", " ", " "]
                        print("\t".join(output))
                else: # HGNC ID not found in uniprot file
                    output = [" ", Candidate_Gene, " ", " ", " ", " ", " ", " "]
                    print("\t".join(output))
        
        print(" ")

    logging.info("All done, completed succesfully!")

    return

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
---------------------------------------------------------------------------------------------------
Program: Parses the Candidate Results File produced by `ExtractCandidates.py` script. Also parses
         the Uniprot Primary Accession File, Canonical Transcripts file and the Orthology data file 
         (retireved from Alliance of Genome Resources, AGR). For each candidate human gene, the 
         script finds a model organism ortholog(s), if available, with additional annotations.

!!!NOTE!!!: This script is written to find Yeast orthologs (TaxID:559292) for human genes.
            If you are interested in finding an ortholog for a different model organism, then you
            will have modify this script (including `1_Uniprot_parses.py` script). These are very 
            minor modifications (particularly TaxID and GeneID lines).
---------------------------------------------------------------------------------------------------
- The output consists of a header followed by data in tab-delimited (.tsv) format. The lines starting
with '#' indicate the name of the confidence level. 
- Following information is provided for each candidate gene:
    # - Name of the Human candidate gene (including Synonyms)
    # - HGNC ID
    # - Whether the candidate gene is already a known candidate for the phenotype interested in
    # - Human candidate gene function
    # - Model organism ortholog (if available)
    # - Model organism Gene ID (as per AGR)
    # - Orthology match scores
    # - Model organsim Gene function
---------------------------------------------------------------------------------------------------

Arguments [defaults] -> Can be abbreviated to shortest unambiguous prefixes
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)


    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inuniprot', metavar = "Input File", dest = "inuniprot", help = 'Uniprot Primary Accession File generated by the uniprot parser', required = True)
    required.add_argument('--incanonical', metavar = "Input File", dest = "incanonical", help = 'Canonical Transcripts file (.gz or tsv)', required = True)
    required.add_argument('--incandidateresult', metavar = "Input File", dest = "incandidateresult", help = 'File containing list of candidates (different confidence levels) extracted by ExtractCandidates.py', required = True)
    required.add_argument('--inortholog', metavar = "Input File", dest = "inortholog", help = 'Orthology data file retireved from Alliance of Genome Resources', required = True)

    args = file_parser.parse_args()
    FindOrthologs(args)


if __name__ == "__main__":
    # Logging to Standard Error
    Log_Format = "%(levelname)s - %(asctime)s - %(message)s \n"
    logging.basicConfig(stream = sys.stderr, format  = Log_Format, level = logging.DEBUG)
    main()
