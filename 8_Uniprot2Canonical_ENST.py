#!/usr/bin/python

import argparse, sys
import re

###Function for creating ranscripts_Gene dictionary###
# Takes tab-seperated canonical transcripts file as INPUT
# Creates a dictionary using 2 columns: ENST and Gene
# Key -> ENST; Value - Gene
# Returns the dictionary
def Transcripts_Gene(args):

    Transcripts_Gene_dict = {} # Initializing an empty dictionary

    CanonicalTranscripts_File = open(args.inTranscripts)

    CanonicalTranscripts_File.readline() # Skip header

    # Parsing the Uniprot Primary Accession file
    for line in CanonicalTranscripts_File:
        line = line.rstrip('\n')
        CanonicalTranscripts_fields = line.split('\t')

        # Key -> ENST
        # Value -> Gene
        (ENST_key, Gene) = (CanonicalTranscripts_fields[0], CanonicalTranscripts_fields[1])
        Transcripts_Gene_dict[ENST_key] = Gene

    return Transcripts_Gene_dict

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

    # Calling the dictionary
    Transcripts_Gene_dict = Transcripts_Gene(args)

    UniprotPrimAC_file = open(args.inPrimAC)

    UniprotPrimAC_file.readline() # Skip header

    Uniprot_ENST_dict = {} # Initializing an empty dictionary

    # Compiling regular expressions###

    # Eliminating Mouse transcripts
    re_ENSMUST = re.compile('^ENSMUST')

    # Counter for canonical transcripts
    canonical_transcripts_count = 0

    # Counter for accessions with no canonical human ENST
    no_CanonicalHumanENST = 0

    # Counter for accessions with single canonical human ENST
    single_CanonicalHumanENST = 0

    # Counter for accessions with multiple canonical human ENSTs
    multiple_CanonicalHumanENST = 0

    # Parsing the Uniprot Primary Accession file
    for line in UniprotPrimAC_file:
        line = line.rstrip('\n')
        UniprotPrimAC_fields = line.split('\t')

        # ENST column - UniprotPrimAC_fields[2]
        # This is a single string containing comma-seperated ENSTs
        # So we split it into a list that can be accessed later
        UniProt_ENSTs = UniprotPrimAC_fields[2].split(',')

        # Initializing an empty ENST list
        human_ENSTs = []
        canonical_human_ENSTs = []

        # Eliminating Mouse transcripts
        for UniProt_ENST in UniProt_ENSTs:
            if not re_ENSMUST.match(UniProt_ENST):
                human_ENSTs.append(UniProt_ENST)

                for ENST in human_ENSTs:
                    if ENST in Transcripts_Gene_dict.keys():
                        canonical_human_ENSTs.append(ENST)
                        canonical_transcripts_count += 1

                # Key -> Uniprot Primary accession
                # Value -> Canonical_human_ENST
                (UniprotPrimAC_key, canonical_ENST) = (UniprotPrimAC_fields[0], canonical_human_ENSTs)

                # After eliminating mouse and non-canonical transcripts, some values can be empty
                # So, we keep a count of these accessions
                if len(canonical_human_ENST) == 0:
                    no_CanonicalHumanENST += 1
                elif len(canonical_human_ENST) == 1:
                    Uniprot_ENST_dict[UniprotPrimAC_key] = canonical_ENST
                    single_CanonicalHumanENST += 1
                elif len(canonical_human_ENST) > 1:
                    multiple_CanonicalHumanENST += 1

    print("\nTotal no. of Canonical transcripts:", canonical_transcripts_count, file = sys.stderr)
    print("\nNo. of UniProt primary accessions without canonical human transcripts:", no_CanonicalHumanENST, file = sys.stderr)
    print("\nNo. of UniProt primary accessions with single canonical human transcript:", single_CanonicalHumanENST, file = sys.stderr)
    print("\nNo. of UniProt primary accessions with multiple canonical human transcripts:", multiple_CanonicalHumanENST, file = sys.stderr)

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
 -> Corresponding ENST
---------------------------------------------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inPrimAC', metavar = "Input File", dest = "inPrimAC", help = 'Uniprot Primary Accession File generated by the uniprot parser', required = True)
    required.add_argument('--inTranscripts', metavar = "Input File", dest = "inTranscripts", help = 'Canonical Transcripts file', required = True)

    file_parser.set_defaults(func=Uniprot2ENST)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
