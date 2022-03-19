#!/usr/bin/python

import sys, argparse
import re
import logging
import time

###########################################################

# Parses the output files produced by the interaction_parser.py
# Required columns are: Protein A UniProt PrimAC, Protein B UniProt PrimAC,
# Interaction Detection Method, Pubmed Identifier and Interaction type
#
# Processes it by filtering based on Interaction Detection Method and Interaction type
# Final filtering - Each interaction has at least 2 experiments
# Returns a list with 5 items (in each sublist):
# - Protein A UniProt PrimAC,
# - Protein B UniProt PrimAC,
# - Publication_Count,
# - PMID(s) and
# - Experiment_count
def UniProtInteractome(inCuratedFile):

    PPI_PMID_dict = {} # Dictionary for PPIs
    PPI_IntDetMethod_dict = {} # Dictionary for experiments and filtering based on IntDetMethod

    # List of User input curated interaction files
    curatedFiles = inCuratedFile

    # there can be multiple files
    for file in curatedFiles:

        logging.info("Processing data from File: %s" % file)

        curatedIntFile = open(file)

        # Parsing the file
        for line in curatedIntFile:

            line = line.rstrip('\n')

            curatedPPI_fields = line.split('\t')

            # Filtering out Interactions based on Interaction Detection Method
            IntDetMethod = curatedPPI_fields[2]
            # MI:0096 -> pull down
            # MI:0254 -> genetic interference
            # MI:0686 -> unspecified method

            ###Filtering out Interactions based on Interaction Type
            IntType = curatedPPI_fields[4].rstrip('\n')
            # MI:0407 -> direct interaction
            # MI:0915 -> physical association

            if IntDetMethod not in ['MI:0096', 'MI:0254', 'MI:0686'] and IntType in ['MI:0407', 'MI:0915']:
                # curatedPPI_fields[0] -> Protein_A_UniprotPrimAC
                # curatedPPI_fields[1] -> Protein_B_UniprotPrimAC
                Interactors = curatedPPI_fields[0] + '_' + curatedPPI_fields[1]

                # Key -> UniProt PrimAC of Protein A & B joined together by an '_'
                # Value -> Pubmed Identifier (PMID) - curatedPPI_fields[3]
                (Int_key, PMIDs) = (Interactors, curatedPPI_fields[3])

                # Check if the Key exists in PPI_PMID_dict
                # If yes, then store the values (PMIDs) as a list
                if PPI_PMID_dict.get(Int_key, False):
                    if PMIDs not in PPI_PMID_dict[Int_key]: # Avoiding duplicate PMIDs
                        PPI_PMID_dict[Int_key].append(PMIDs)
                else:
                    PPI_PMID_dict[Int_key] = [PMIDs]

                # Key -> UniProt PrimAC of Protein A & B joined together by an '_'
                # Value -> Interaction Detection Method - curatedPPI_fields[2]
                (Int_key, IntDetMeth) = (Interactors, curatedPPI_fields[2])

                if PPI_IntDetMethod_dict.get(Int_key, False):
                    PPI_IntDetMethod_dict[Int_key].append(IntDetMeth)
                else:
                    PPI_IntDetMethod_dict[Int_key] = [IntDetMeth]

        # Closing the file
        curatedIntFile.close()

    # Initializing output list
    Uniprot_Interactome_list = []

    logging.info("Building High-Quality Human Interactome...")

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
            if Exp_count >= 2:
                interaction_out_line = [Protein_A, Protein_B, str(PMID_count), Pubmed_Identifier, str(Exp_count)]
                Uniprot_Interactome_list.append(interaction_out_line)

    return Uniprot_Interactome_list

###########################################################

# Parses tab-seperated canonical transcripts file
# Required columns are: 'ENSG' and 'GENE' (can be in any order,
# but they MUST exist)
#
# Returns a dictionary:
# Key -> ENSG
# Value -> Gene
def ENSG_Gene(inCanonicalFile):

    ENSG_Gene_dict = {} # Initializing an empty dictionary

    Canonical_File = open(inCanonicalFile)

    Canonical_header_line = Canonical_File.readline() # Grabbing the header line

    Canonical_header_fields = Canonical_header_line.split('\t')

    # Check the column headers and grab indexes of our columns of interest
    (ENSG_col, Gene_col) = (-1,-1)

    logging.info("Processing data from Canonical Transcripts File: %s" % inCanonicalFile)

    for i in range(len(Canonical_header_fields)):
        if Canonical_header_fields[i] == 'ENSG':
            ENSG_col = i
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

        ENSG_Gene_dict[CanonicalTranscripts_fields[ENSG_col]] = CanonicalTranscripts_fields[Gene_col]

    # Closing the file
    Canonical_File.close()

    return ENSG_Gene_dict

###########################################################

# Parses the UniProt Primary Accession file produced by uniprot_parser.py
# Required columns are: 'Primary_AC' and 'ENSGs' (can be in any order,
# but they MUST exist)
#
# Parses the dictionary ENSG_Gene_dict
# returned by the function ENSG_Gene
# Maps UniProt Primary Accession to ENSG
# Returns a dictionary:
# - Key: UniProt Primary Accession
# - Value: Corresponding ENSG

def Uniprot_ENSG(inPrimAC, ENSG_Gene_dict):

    # Initializing the dictionary
    Uniprot_ENSG_dict = {}

    UniprotPrimAC_File = open(inPrimAC)

    # Grabbing the header line
    UniprotPrimAC_header = UniprotPrimAC_File.readline()

    UniprotPrimAC_header = UniprotPrimAC_header.rstrip('\n')

    UniprotPrimAC_header_fields = UniprotPrimAC_header.split('\t')

    # Check the column header and grab indexes of our columns of interest
    (UniProt_PrimAC_col, ENSG_col) = (-1, -1)

    logging.info("Processing data from Uniprot Primary Accession File: %s" % inPrimAC)

    for i in range(len(UniprotPrimAC_header_fields)):
        if UniprotPrimAC_header_fields[i] == 'Primary_AC':
            UniProt_PrimAC_col = i
        elif UniprotPrimAC_header_fields[i] == 'ENSGs':
            ENSG_col = i

    if not UniProt_PrimAC_col >= 0:
        sys.exit("Missing required column title 'Primary_AC' in the file: %s \n" % inPrimAC)
    elif not ENSG_col >= 0:
        sys.exit("Missing required column title 'ENSG' in the file: %s \n" % inPrimAC)
    # else grabbed the required column indexes -> PROCEED

    # Compiling regular expressions###

    # Eliminating Mouse ENSGs
    re_ENSMUST = re.compile('^ENSMUSG')

    #Counter for Human Uniprot Primary Accession
    Count_HumanUniprotPrimAC = 0

    # Counter for canonical ENSGs
    canonical_ENSG_count = 0

    # Counter for accessions with no canonical human ENSG
    no_CanonicalHumanENSG = 0

    # Counter for accessions with single canonical human ENSG
    single_CanonicalHumanENSG = 0

    # Counter for accessions with multiple canonical human ENSGs
    multiple_CanonicalHumanENSG = 0

    # Parsing the Uniprot Primary Accession file
    for line in UniprotPrimAC_File:
        line = line.rstrip('\n')
        UniprotPrimAC_fields = line.split('\t')

        # ENSG column  - This is a single string containing comma-seperated ENSGs
        # So we split it into a list that can be accessed later
        UniProt_ENSGs = UniprotPrimAC_fields[ENSG_col].split(',')

        # Initializing empty lists
        human_ENSGs = []
        canonical_human_ENSGs = []

        # Eliminating Mouse ENSGs
        for UniProt_ENSG in UniProt_ENSGs:
            if not re_ENSMUST.match(UniProt_ENSG):
                human_ENSGs.append(UniProt_ENSG)
        if not human_ENSGs:
            continue
        Count_HumanUniprotPrimAC += 1

        for ENSG in human_ENSGs:
            if ENSG in ENSG_Gene_dict.keys():
                canonical_human_ENSGs.append(ENSG)
                canonical_ENSG_count += 1

        # Key -> Uniprot Primary accession
        # Value -> Canonical_human_ENSG

        # After eliminating mouse and non-canonical ENSGs, some values can be empty
        # So, we keep a count of these accessions
        if len(canonical_human_ENSGs) == 0:
            no_CanonicalHumanENSG += 1
        elif len(canonical_human_ENSGs) == 1:
            Uniprot_ENSG_dict[UniprotPrimAC_fields[UniProt_PrimAC_col]] = ''.join(canonical_human_ENSGs)
            single_CanonicalHumanENSG += 1
        elif len(canonical_human_ENSGs) > 1:
            multiple_CanonicalHumanENSG += 1

    logging.debug("Total no. of Human UniProt Primary Accessions: %d " % Count_HumanUniprotPrimAC)
    logging.debug("Total no. of ENSGs in the Canonical Transcripts file: %d " % len(ENSG_Gene_dict.keys()))
    logging.debug("Total no. of ENSGs in the UniProt Primary Accession file: %d " % canonical_ENSG_count)
    logging.debug("No. of UniProt primary accessions without canonical human ENSG: %d " % no_CanonicalHumanENSG)
    logging.debug("No. of UniProt primary accessions with single canonical human ENSG: %d " % single_CanonicalHumanENSG)
    logging.debug("No. of UniProt primary accessions with multiple canonical human ENSGs: %d " % multiple_CanonicalHumanENSG)

    # Closing the file
    UniprotPrimAC_File.close()

    return Uniprot_ENSG_dict

###########################################################

# Parses the list - [Interactome_list] returned by the function: UniProtInteractome and
# the dictionary - {Uniprot_ENSG_dict} returned by the function: Uniprot_ENSG
#
# Maps the UniProt Primary Accessions to ENSG using the dictionary: {ENSG_Gene_dict}
# Prints to STDOUT in .tsv format
# Output consists of 5 columns:
# - ENSG of Protein A
# - ENSG of Protein B
# - Number of Publications associated with the interaction of the above 2 proteins
# - PMID (or a comma seperated list of PMIDs)
# - Count of Experiments for each interaction
def Interactome_Uniprot2ENSG(args):

    # Calling the functions
    Uniprot_Interactome_list = UniProtInteractome(args.inCuratedFile)
    ENSG_Gene_dict = ENSG_Gene(args.inCanonicalFile)
    Uniprot_ENSG_dict = Uniprot_ENSG(args.inPrimAC, ENSG_Gene_dict)

    # Counter for UniProt Primary Accessions of proteins not mapping to ENSG
    lost_Interaction = 0

    logging.info("Mapping UniProt Primary Accessions to ENSG")

    for data in Uniprot_Interactome_list:

        if data[0] in Uniprot_ENSG_dict.keys() and data[1] in Uniprot_ENSG_dict.keys():
            ENSG_Interactome_out = (Uniprot_ENSG_dict.get(data[0]), Uniprot_ENSG_dict.get(data[1]), data[2], data[3], data[4])
            print('\t'.join(ENSG_Interactome_out))
        else:
            lost_Interaction += 1


    logging.debug("Total no. of Interactions lost: %d " % lost_Interaction)
    logging.info("Done 🎉")

    return

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
-------------------------------------------------------------------------------------------------------------
Program: Parses the output file(s) produced by the interaction_parser.py to produce a high-quality interactome,
         maps the UniProt Primary Accession of interacting proteins to ENSG in the Canonical transcripts file
         and prints to STDOUT
-------------------------------------------------------------------------------------------------------------

The output (High-quality Human Interactome) consists of five columns in .tsv format:
  -> ENSG of Protein A
  -> ENSG of Protein B
  -> Number of Publications associated with the interaction of the above 2 proteins
  -> PMID (or a comma seperated list of PMIDs)
  -> Count of Experiments for each interaction
--------------------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)


    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inCuratedFile', metavar = "Input File", dest = "inCuratedFile", nargs = '+', help = 'Output files produced by interaction_parser.py', required = True)
    required.add_argument('--inPrimAC', metavar = "Input File", dest = "inPrimAC", help = 'Uniprot Primary Accession File generated by the uniprot parser', required = True)
    required.add_argument('--inCanonicalFile', metavar = "Input File", dest = "inCanonicalFile", help = 'Canonical Transcripts file', required = True)

    args = file_parser.parse_args()
    Interactome_Uniprot2ENSG(args)


if __name__ == "__main__":
    # Logging
    Log_Format = "%(levelname)s - %(asctime)s - %(message)s \n"
    logging.basicConfig(stream = sys.stderr, format  = Log_Format, level = logging.DEBUG)
    main()
