#!/usr/bin/python

import sys, argparse
import logging

###########################################################

# Parses tab-seperated canonical transcripts file
# Required columns are: 'ENSG' and 'GENE' (can be in any order,
# but they MUST exist)
# Returns a dictionary:
# Key -> ENSG; Value -> Gene
def ENSG_Gene(inCanonicalFile):

    logging.info("Starting to run...")

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
        sys.exit("Missing required column title: 'ENSG' in the file: %s \n" % inCanonicalFile)
    elif not Gene_col >= 0:
        sys.exit("Missing required column title: 'GENE' in the file: %s \n" % inCanonicalFile)
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

###########################################################

# Parses the UniProt Primary Accession file produced by uniprot_parser.py
# Required columns are: 'Primary_AC' and 'ENSG' (can be in any order,
# but they MUST exist)
#
# Parses the dictionary returned by the function ENSG_Gene
# Maps UniProt Primary Accession to ENSG
# Prints to STDOUT in .tsv format
# Output consists of 2 columns in .tsv format:
# - Uniprot Primary Accession
# - Corresponding ENSG
def Uniprot2ENSG(args):

    # Calling the function ENSG_Gene
    ENSG_Gene_dict = ENSG_Gene(args.inCanonicalFile)

    Uniprot_File = open(args.inUniProt)

    # Grabbing the header line
    Uniprot_header = Uniprot_File.readline()

    Uniprot_header = Uniprot_header.rstrip('\n')

    Uniprot_header_fields = Uniprot_header.split('\t')

    # Check the column headers and grab indexes of our columns of interest
    (UniProt_PrimAC_index, ENSG_index) = (-1, -1)

    for i in range(len(Uniprot_header_fields)):
        if Uniprot_header_fields[i] == 'Primary_AC':
            UniProt_PrimAC_index = i
        elif Uniprot_header_fields[i] == 'ENSGs':
            ENSG_index = i

    if not UniProt_PrimAC_index >= 0:
        logging.error("At Step 5.2_addInteractome - Missing required column title 'Primary_AC' in the file: %s \n" % inUniProt)
        sys.exit()
    elif not ENSG_index >= 0:
        logging.error("At Step 5.2_addInteractome - Missing required column title 'ENSG' in the file: %s \n" % inUniProt)
        sys.exit()
    # else grabbed the required column indexes -> PROCEED


    # Data lines
    for line in Uniprot_File:
        line = line.rstrip('\n')
        Uniprot_fields = line.split('\t')

        # ENSG column  - This is a single string containing comma-seperated ENSGs
        # So we split it into a list that can be accessed later
        UniProt_ENSGs = Uniprot_fields[ENSG_index].split(',')
        
        canonical_human_ENSGs = []

        # If ENSG is in the canonical transcripts file
        # Append it to canonical_human_ENSGs
        for ENSG in UniProt_ENSGs:
            if ENSG in ENSG_Gene_dict.keys():
                canonical_human_ENSGs.append(ENSG)

        # Keeping the count of protein with single ENSGs
        if len(canonical_human_ENSGs) == 1:
            print(Uniprot_fields[UniProt_PrimAC_index], "\t", ''.join(canonical_human_ENSGs))

    Uniprot_File.close()

    logging.info("All done, completed successfully!")

    return

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
---------------------------------------------------------------------------------------------------------------------------------
Program: Parses the Uniprot file (produced by 1_Uniprot_parser.py) and the canonical transcripts file, maps the Uniprot Primary 
         Accessions to ENSG and prints to STDOUT 
---------------------------------------------------------------------------------------------------------------------------------
The output consists of 2 columns in .tsv format:
 -> Uniprot Primary Accession
 -> Corresponding ENSG
---------------------------------------------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inUniProt', metavar = "Input File", dest = "inUniprot", help = 'Uniprot output File generated by the UniProt_parser.py', required = True)
    required.add_argument('--inCanonicalFile', metavar = "Input File", dest = "inCanonicalFile", help = 'Canonical Transcripts file', required = True)

    args = file_parser.parse_args()
    Uniprot2ENSG(args)

if __name__ == "__main__":
    # Logging to Standard Error
    Log_Format = "%(levelname)s - %(asctime)s - %(message)s \n"
    logging.basicConfig(stream = sys.stderr, format  = Log_Format, level = logging.DEBUG)
    main()
