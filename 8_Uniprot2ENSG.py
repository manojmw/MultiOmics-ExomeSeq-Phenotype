#!/usr/bin/python

import argparse, sys
import re

###Function for creating transcripts_Gene dictionary###
# Takes tab-seperated canonical transcripts file as INPUT
# Creates a dictionary using 2 columns: ENSG and Gene
# Key -> ENSG; Value - Gene
# Returns the dictionary
def ENSG_Gene(inCanonicalFile):

    ENSG_Gene_dict = {} # Initializing an empty dictionary

    Canonical_File = open(inCanonicalFile)

    Canonical_header_line = Canonical_File.readline() # Grabbing the header line

    Canonical_header_fields = Canonical_header_line.split('\t')

    # Check the column headers and grab indexes of our columns of interest
    (ENSG_col, Gene_col) = (-1,-1)

    for header in Canonical_header_fields:
        if header == 'ENSG':
            ENSG_col = Canonical_header_fields.index(header)
        elif header == 'GENE':
            Gene_col = Canonical_header_fields.index(header)

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
        (ENSG_key, Gene) = (CanonicalTranscripts_fields[ENSG_col], CanonicalTranscripts_fields[Gene_col])
        ENSG_Gene_dict[ENSG_key] = Gene

    return ENSG_Gene_dict

###Function for mapping Uniprot Primary Accession to canonical transcripts###
# Calls the above dictionary - Transcripts_Gene_dict
# Takes Uniprot Primary Accession file produced by uniprot_parser as INPUT
# Creates another dictionary - Uniprot_ENST_dict using 2 columns
# i.e. Uniprot Primary accession and ENST(s); # Key -> Uniprot Primary accession; Value - ENST(s)
# Using these dictionaries, checks whether ENST(s) of each accession from Uniprot Primary Accession file
# is present in the canonical transcripts file
# If yes, checks whether this accession has single canonical_ENST asssociated
# If yes, print the Accession and ENST to STDOUT in .tsv format
# If multiple canonical_ENSTs or no canonical_ENST for a given UniProt accession
# Then, do not print and keep a count seperately
def Uniprot2ENST(args):

    # Calling the function ENSG_Gene
    Transcripts_Gene_dict = ENSG_Gene(args.inCanonicalFile)

    # INPUT -> Uniprot Primary_AC File
    UniprotPrimAC_File = open(args.inPrimAC)

    UniprotPrimACFile_header_line = UniprotPrimAC_File.readline() # Grabbing the header line

    UniprotPrimACFile_header_line = UniprotPrimACFile_header_line.strip('\n')

    UniprotPrimACFile_header_fields = UniprotPrimACFile_header_line.split('\t')

    # Check the column headers and grab indexes of our columns of interest
    (UniProt_PrimAC_col, ENSG_col) = (-1,-1)

    for header in UniprotPrimACFile_header_fields:
        if header == 'Primary_AC':
            UniProt_PrimAC_col = UniprotPrimACFile_header_fields.index(header)
        elif header == 'ENSG':
            ENSG_col = UniprotPrimACFile_header_fields.index(header)

    if not UniProt_PrimAC_col >= 0:
        sys.exit("Missing required column title: 'Primary_AC' \n")
    elif not ENSG_col >= 0:
        sys.exit("Missing required column title: 'ENSG' \n")
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
            if ENSG in Transcripts_Gene_dict.keys():
                canonical_human_ENSGs.append(ENSG)
                canonical_ENSG_count += 1

        # Key -> Uniprot Primary accession
        # Value -> Canonical_human_ENSG

        # After eliminating mouse and non-canonical ENSGs, some values can be empty
        # So, we keep a count of these accessions
        if len(canonical_human_ENSGs) == 0:
            no_CanonicalHumanENSG += 1
        elif len(canonical_human_ENSGs) == 1:
            print(UniprotPrimAC_fields[UniProt_PrimAC_col], '\t', ''.join(canonical_human_ENSGs))
            single_CanonicalHumanENSG += 1
        elif len(canonical_human_ENSGs) > 1:
            multiple_CanonicalHumanENSG += 1

    print("\nTotal no. of Human UniProt Primary Accessions:", Count_HumanUniprotPrimAC, file = sys.stderr)
    print("\nTotal no. of ENSGs in the Canonical Transcripts file:", len(Transcripts_Gene_dict.keys()), file = sys.stderr)
    print("\nTotal no. of ENSGs in the UniProt Primary Accession file:", canonical_ENSG_count, file = sys.stderr)
    print("\nNo. of UniProt primary accessions without canonical human ENSG:", no_CanonicalHumanENSG, file = sys.stderr)
    print("\nNo. of UniProt primary accessions with single canonical human ENSG:", single_CanonicalHumanENSG, file = sys.stderr)
    print("\nNo. of UniProt primary accessions with multiple canonical human ENSGs:", multiple_CanonicalHumanENSG, file = sys.stderr)

    return

####Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
---------------------------------------------------------------------------------------------------------------------------------
Program: Parses the Uniprot Primary Accession file produced by uniprot_parser, maps to canonical transcripts and prints to STDOUT
---------------------------------------------------------------------------------------------------------------------------------
The output consists of 2 columns in .tsv format:
 -> Uniprot Primary Accession
 -> Corresponding ENSG
---------------------------------------------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inPrimAC', metavar = "Input File", dest = "inPrimAC", help = 'Uniprot Primary Accession File generated by the uniprot parser', required = True)
    required.add_argument('--inCanonicalFile', metavar = "Input File", dest = "inCanonicalFile", help = 'Canonical Transcripts file', required = True)

    args = file_parser.parse_args()
    Uniprot2ENST(args)

if __name__ == "__main__":
    main()
