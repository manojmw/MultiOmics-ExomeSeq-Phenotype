#!/usr/bin/python

# manojmw
# 27 Jan, 2022

import sys, argparse
import logging
import gzip

###########################################################

# Parses the tab-seperated output files
# produced by the Interaction_parser.py (No headers)
#
# Required columns are:
# - Protein A UniProt PrimAC (column-1)
# - Protein B UniProt PrimAC (column-2)
# - Interaction Detection Method (column-3)
# - Pubmed Identifier (column-4)
# - Interaction type (column-5)
#
# Processes it by filtering based on Interaction Detection Method and Interaction type
# Final filtering - Each interaction has at least 2 experiments
# at least one of the experiments proved by any binary interaction method
#
# Returns a list with 5 items (in each sublist):
# - Protein A UniProt PrimAC,
# - Protein B UniProt PrimAC,
# - Publication_Count,
# - PMID(s) and
# - Experiment_count
def UniProtInteractome(inExpFile):

    logging.info("Starting to run...")

    # Dictionary for storing Interacting
    # Proteins and Pubmed Identifier
    PPI_PMID_dict = {}

    # Dictionary for storing Interacting
    # Proteins and Interaction Detection Method
    PPI_IntDetMethod_dict = {}

    # List of User input PPI interaction experiments file(s)
    PPIExpFiles = inExpFile

    # there can be multiple files
    for file in PPIExpFiles:

        PPIExpFile = open(file)

        # Data lines
        for line in PPIExpFile:

            line = line.rstrip('\n')

            ExpPPI_fields = line.split('\t')

            # Filtering out Interactions based on Interaction Detection Method
            IntDetMethod = ExpPPI_fields[2]
            # MI:0096 -> pull down
            # MI:0254 -> genetic interference
            # MI:0686 -> unspecified method

            ###Filtering out Interactions based on Interaction Type
            IntType = ExpPPI_fields[4].rstrip('\n')
            # MI:0407 -> direct interaction
            # MI:0915 -> physical association

            if IntDetMethod not in ['MI:0096', 'MI:0254', 'MI:0686'] and IntType in ['MI:0407', 'MI:0915']:
                # ExpPPI_fields[0] -> Protein_A_UniprotPrimAC
                # ExpPPI_fields[1] -> Protein_B_UniprotPrimAC
                Interactors = ExpPPI_fields[0] + '_' + ExpPPI_fields[1]

                # Key -> UniProt PrimAC of Protein A & B joined together by an '_'
                # Value -> Pubmed Identifier (PMID) - ExpPPI_fields[3]
                (Int_key, PMIDs) = (Interactors, ExpPPI_fields[3])

                # Check if the Key exists in PPI_PMID_dict
                # If yes, then store the values (PMIDs) as a list
                if PPI_PMID_dict.get(Int_key, False):
                    # Avoiding duplicate PMIDs
                    if PMIDs not in PPI_PMID_dict[Int_key]:
                        PPI_PMID_dict[Int_key].append(PMIDs)
                else:
                    PPI_PMID_dict[Int_key] = [PMIDs]

                # Key -> UniProt PrimAC of Protein A & B joined together by an '_'
                # Value -> Interaction Detection Method - ExpPPI_fields[2]
                (Int_key, IntDetMeth) = (Interactors, ExpPPI_fields[2])

                if PPI_IntDetMethod_dict.get(Int_key, False):
                    PPI_IntDetMethod_dict[Int_key].append(IntDetMeth)
                else:
                    PPI_IntDetMethod_dict[Int_key] = [IntDetMeth]

        # Closing the file
        PPIExpFile.close()

    # Initializing output list
    Uniprot_Interactome_list = []

    # Processing the dictionaries and returning a list
    # Since both the dictionaries have the same keys,
    # We can use the keys from PPI_PMID_dict to iterate through PPI_IntDetMethod_dict
    for Int_key in PPI_PMID_dict:

        # Checking if at least one of the experiments for a given interaction has been proved by any binary interaction method
        # i.e. Other than Affintity Chromatography Technology (ACT) - 'MI:0004'
        # Used to eliminate PPIs that have been proved ONLY by using ACT
        if any(value != 'MI:0004' for value in PPI_IntDetMethod_dict[Int_key]):
            Proteins = Int_key.split('_')
            Protein_A = Proteins[0]
            Protein_B = Proteins[1]

            Pubmed_Identifier = ', '.join(PPI_PMID_dict[Int_key])
            PMID_count = len(PPI_PMID_dict[Int_key])
            Exp_count = len(PPI_IntDetMethod_dict[Int_key])

            # Final Quality Control
            # Each interaction has at least 2 experiments
            if (Exp_count >= 2):
                interaction_out_line = [Protein_A, Protein_B, str(PMID_count), Pubmed_Identifier, str(Exp_count)]
                Uniprot_Interactome_list.append(interaction_out_line)

    return Uniprot_Interactome_list

###########################################################

# Parses tab-seperated canonical transcripts file
# Required columns are: 'ENSG' and 'GENE' (can be in any order,
# but they MUST exist)
#
# Returns a dictionary:
# - Key -> ENSG
# - Value -> Gene
def ENSG_Gene(inCanonicalFile):

    # Dictionary for storing ENSG
    # and Gene data
    ENSG_Gene_dict = {}

     # Opening canonical transcript file (gzip or non-gzip)
    try:
        if inCanonicalFile.endswith('.gz'):
            Canonical_File = gzip.open(inCanonicalFile, 'rt')
        else:
            Canonical_File = open(inCanonicalFile)
    except IOError:
        sys.exit("Error: Failed to read the Canonical transcript file: %s" % inCanonicalFile)

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
        sys.exit("Error: Missing required column title 'ENSG' in the file: %s \n" % inCanonicalFile)
    elif not Gene_index >= 0:
        sys.exit("Error: Missing required column title 'GENE' in the file: %s \n" % inCanonicalFile)
    # else grabbed the required column indexes -> PROCEED

    # Data lines
    for line in Canonical_File:
        line = line.rstrip('\n')
        CanonicalTranscripts_fields = line.split('\t')

        # Key -> ENSG
        # Value -> Gene

        ENSG_Gene_dict[CanonicalTranscripts_fields[ENSG_index]] = CanonicalTranscripts_fields[Gene_index]

    # Closing the file
    Canonical_File.close()

    return ENSG_Gene_dict

###########################################################

# Parses the tab-seperated UniProt Primary Accession file
# produced by Uniprot_parser.py
#
# Required columns are: 'Primary_AC' and 'ENSGs' (can be in any order,
# but they MUST exist)
#
# Parses the dictionary ENSG_Gene_dict
# returned by the function ENSG_Gene
# Maps UniProt Primary Accession to ENSG
#
# Returns a dictionary:
# - Key: UniProt Primary Accession
# - Value: Corresponding ENSG

def Uniprot_ENSG(inUniProt, ENSG_Gene_dict):

    # Dictionary for storing UniProt Primary
    # Accession and ENSG data
    Uniprot_ENSG_dict = {}

    Uniprot_File = open(inUniProt)

    # Grabbing the header line
    Uniprot_header = Uniprot_File.readline()

    Uniprot_header = Uniprot_header.rstrip('\n')

    Uniprot_header_fields = Uniprot_header.split('\t')

    # Check the column header and grab indexes of our columns of interest
    (UniProt_PrimAC_index, ENSG_index) = (-1, -1)

    for i in range(len(Uniprot_header_fields)):
        if Uniprot_header_fields[i] == 'Primary_AC':
            UniProt_PrimAC_index = i
        elif Uniprot_header_fields[i] == 'ENSGs':
            ENSG_index = i

    # Sanity check
    if not UniProt_PrimAC_index >= 0:
        sys.exit("Error: Missing required column title 'Primary_AC' in the file: %s \n" % inUniProt)
    elif not ENSG_index >= 0:
        sys.exit("Error: Missing required column title 'ENSG' in the file: %s \n" % inUniProt)
    # else grabbed the required column indices -> PROCEED

    # Counter for canonical ENSGs
    canonical_ENSG_count = 0

    # Counter for accessions with no canonical human ENSG
    no_CanonicalHumanENSG = 0

    # Counter for accessions with single canonical human ENSG
    single_CanonicalHumanENSG = 0

    # Counter for accessions with multiple canonical human ENSGs
    multiple_CanonicalHumanENSG = 0

    # Data lines
    for line in Uniprot_File:
        line = line.rstrip('\n')
        Uniprot_fields = line.split('\t')

        # ENSG column  - This is a single string containing comma-seperated ENSGs
        # So we split it into a list that can be accessed later
        UniProt_ENSGs = Uniprot_fields[ENSG_index].split(',')

        # List to store Human ENSG(s)
        # found in the canonical transcripts file
        canonical_human_ENSGs = []

        for ENSG in UniProt_ENSGs:
            if ENSG in ENSG_Gene_dict.keys():
                canonical_human_ENSGs.append(ENSG)
                canonical_ENSG_count += 1

        # Key -> Uniprot Primary accession
        # Value -> Canonical_human_ENSG

        if len(canonical_human_ENSGs) == 0:
            no_CanonicalHumanENSG += 1
        elif len(canonical_human_ENSGs) == 1:
            Uniprot_ENSG_dict[Uniprot_fields[UniProt_PrimAC_index]] = ''.join(canonical_human_ENSGs)
            single_CanonicalHumanENSG += 1
        # Also counting accessions with multiple ENSGs
        elif len(canonical_human_ENSGs) > 1:
            multiple_CanonicalHumanENSG += 1

    # logging.debug("Total no. of Human UniProt Primary Accessions: %d " % Count_HumanUniprotPrimAC)
    # logging.debug("Total no. of ENSGs in the Canonical Transcripts file: %d " % len(ENSG_Gene_dict.keys()))
    # logging.debug("Total no. of ENSGs in the UniProt Primary Accession file: %d " % canonical_ENSG_count)
    # logging.debug("No. of UniProt primary accessions without canonical human ENSG: %d " % no_CanonicalHumanENSG)
    # logging.debug("No. of UniProt primary accessions with single canonical human ENSG: %d " % single_CanonicalHumanENSG)
    # logging.debug("No. of UniProt primary accessions with multiple canonical human ENSGs: %d " % multiple_CanonicalHumanENSG)

    # Closing the file
    Uniprot_File.close()

    return Uniprot_ENSG_dict

###########################################################

# Parses the Uniprot_Interactome_list returned by 
# the function UniProtInteractome
#
# Required items in each sublist:
# UniProt Primary Accession of Protein A (first item)
# UniProt Primary Accession of Protein B (second item)
#
# Returns 2 dictionaries
# First dictionary contains:
# - key: Protein A; Value: List of interactors
# Second dictionary contains:
# - key: Protein B; Value: List of interactors
def Interactome_dict(Uniprot_Interactome_list):

    # Dictionaries to store interacting proteins
    # In ProtA_dict, key -> Protein A; Value -> Protein B
    # In ProtB_dict, key -> Protein B; Value -> Protein A
    ProtA_dict = {}
    ProtB_dict = {}
        
    # Data lines
    for data in Uniprot_Interactome_list:

        if data[0] != data[1]: # checking self-interaction
            # Check if the Key(ProtA) exists in ProtA_dict
            # If yes, then append the interctor to 
            # the list of values (Interactors)
            if ProtA_dict.get(data[0], False):
                ProtA_dict[data[0]].append(data[1])
            else:
                ProtA_dict[data[0]] = [data[1]]

            # Check if the Key(ProtB) exists in ProtB_dict
            # If yes, then append the interctor to 
            # the list of values (Interactors)
            if ProtB_dict.get(data[1], False):
                ProtB_dict[data[1]].append(data[0])
            else:
                ProtB_dict[data[1]] = [data[0]]    

        # else:
            # NOOP -> The interaction is a self-interaction

    return ProtA_dict, ProtB_dict


# Parses the ProtA_dict & ProtB_dict dictionaries
# returned by the function: Interactome_dict
#
# Checks if a protein is hub/sticky protein
# The criteria for considering a given protein
# as a hub is if it interacts with > 120 proteins
# This number is based on the degree distribution
# of all the proteins in the entire high-quality Interactome
# before eliminating hub/sticky proteins
#
# Returns a dictionary:
# - Key: UniProt Primary Accession of the Hub protein
# - Value: 1
def getHubProteins(ProtA_dict, ProtB_dict):

    # Initializing dictionary to store hub proteins
    HubProteins = {}

    # Checking the number of Interactors for a 
    # given protein
    for protein in ProtA_dict:
        # Get the no. of Interactors
        InteractorsCount_ProtA_dict = len(ProtA_dict[protein])
        # Check if this protein is also
        # present in ProtB_dict
        if protein in ProtB_dict:
            InteractorsCount_ProtB_dict = 0
            # If present, loop through each interactor to 
            # make sure that the Interactors
            # for the current protein in ProtB_dict was not already seen
            # in the Interactor list of ProtA_dict
            for interactor in ProtB_dict[protein]:
                if not interactor in ProtA_dict[protein]:
                    InteractorsCount_ProtB_dict += 1
            Total_InteractorsCount = InteractorsCount_ProtA_dict + InteractorsCount_ProtB_dict
        # if the protien not present in ProtB_dict
        # simply get the Interactors count from ProtA_dict
        else:
            Total_InteractorsCount = InteractorsCount_ProtA_dict

        # If the protein has > 120 Interactors 
        # it is considered a hub/sticky protein
        # append it to the HubProteins list
        if Total_InteractorsCount > 120:
            HubProteins[protein] = 1

    # Checking the interactor count of proteins
    # present only in ProtB_dict
    for protein in ProtB_dict:
        if not protein in ProtA_dict:
            # checking for Hub protein
            if len(ProtB_dict[protein]) > 120:
                HubProteins[protein] = 1

    return HubProteins 

###########################################################

# Parses the list - [Uniprot_Interactome_list]
# returned by the function: UniProtInteractome
# Also parses the dictionary - {Uniprot_ENSG_dict} returned by the function: Uniprot_ENSG
# and {ENSG_Gene_dict} returned by the function: ENSG_Gene
#
# Maps the UniProt Primary Accessions to ENSG
#
# Prints to STDOUT in .tsv format
# Output consists of 2 columns:
# - ENSG of Protein A
# - ENSG of Protein B
def Interactome_Uniprot2ENSG(args):

    # Calling the functions
    Uniprot_Interactome_list = UniProtInteractome(args.inExpFile)
    ENSG_Gene_dict = ENSG_Gene(args.inCanonicalFile)
    Uniprot_ENSG_dict = Uniprot_ENSG(args.inUniProt, ENSG_Gene_dict)
    (ProtA_dict, ProtB_dict) = Interactome_dict(Uniprot_Interactome_list)
    HubProteins = getHubProteins(ProtA_dict, ProtB_dict)

    # Counter for UniProt Primary Accessions of proteins not mapping to ENSG
    # lost_Interaction = 0

    for data in Uniprot_Interactome_list:
        # Eliminating hub/sticky proteins
        if (data[0] in HubProteins) or (data[1] in HubProteins):
            pass 
        else:
            # Get the ENSG for the UniProt Primary Accessions
            if data[0] in Uniprot_ENSG_dict.keys() and data[1] in Uniprot_ENSG_dict.keys():
                # We have eliminated self-interactions using UniProt Accessions
                # But 2 different UniProt Accessions can still have same ENSG ID
                # Eliminate such self-interactions
                if not (Uniprot_ENSG_dict.get(data[0]) == Uniprot_ENSG_dict.get(data[1])):
                    # Sort
                    if Uniprot_ENSG_dict.get(data[0]) < Uniprot_ENSG_dict.get(data[1]):
                        ENSG_Interactome_out = (Uniprot_ENSG_dict.get(data[0]), Uniprot_ENSG_dict.get(data[1]))
                        print('\t'.join(ENSG_Interactome_out))
                    else:
                        ENSG_Interactome_out = (Uniprot_ENSG_dict.get(data[1]), Uniprot_ENSG_dict.get(data[0]))
                        print('\t'.join(ENSG_Interactome_out))
            #   else: self-interaction -> NOOP
        #   else:
                # lost_Interaction += 1

    # logging.debug("Total no. of Interactions lost: %d " % lost_Interaction)
    logging.info("All done, completed succesfully!")

    return

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
--------------------------------------------------------------------------------------------------
Program: Parses the output file(s) produced by the interaction_parser.py to produce a high-quality
         interactome, maps the UniProt Primary Accession of interacting proteins to ENSG using the
         Canonical transcripts file and prints to STDOUT
--------------------------------------------------------------------------------------------------
The output (High-quality Human Interactome) consists of two columns in .tsv format:
  -> ENSG of Protein A
  -> ENSG of Protein B
--------------------------------------------------------------------------------------------------

Arguments [defaults] -> Can be abbreviated to shortest unambiguous prefixes
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)


    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inExpFile', metavar = "Input File", dest = "inExpFile", nargs = '+', help = 'Output files produced by interaction_parser.py', required = True)
    required.add_argument('--inUniprot', metavar = "Input File", dest = "inUniProt", help = 'Uniprot Primary Accession File generated by the uniprot parser', required = True)
    required.add_argument('--inCanonicalFile', metavar = "Input File", dest = "inCanonicalFile", help = 'Canonical Transcripts file (.gz or non .gz)', required = True)

    args = file_parser.parse_args()
    Interactome_Uniprot2ENSG(args)


if __name__ == "__main__":
    # Logging to Standard Error
    Log_Format = "%(levelname)s - %(asctime)s - %(message)s \n"
    logging.basicConfig(stream = sys.stderr, format  = Log_Format, level = logging.DEBUG)
    main()
