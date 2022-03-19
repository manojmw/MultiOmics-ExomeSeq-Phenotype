#!/usr/bin/python

# manojmw
# 14 Jan, 2022

import re
import sys, argparse
import logging
import time

###########################################################

# Parses the tab-seperated UniProt Primary Accession
# file produced by Uniprot_parser.py
#
# Required columns are: 'Primary_AC' and 'ENSG' (can be in any order,
# but they MUST exist)
#
# Returns a dictionary:
# - Key -> UniProt Primary Accession
# - Value -> TaxID
def PrimAC(inPrimAC):

    # Dictionary to store UniProt Primary Accession
    # and Taxonomy Identifier
    UniprotPrimAC_dict = {}

    UniprotPrimAC_File = open(inPrimAC)

    logging.info("Processing data from UniProt Primary Accession File: %s" % inPrimAC)

    # Grabbing the header line
    UniprotPrimAC_header = UniprotPrimAC_File.readline()

    UniprotPrimAC_header = UniprotPrimAC_header.rstrip('\n')

    UniprotPrimAC_header_fields = UniprotPrimAC_header.split('\t')

    # Check the column header and grab indexes of our columns of interest
    (UniProt_PrimAC_index, TaxID_index) = (-1, -1)

    for i in range(len(UniprotPrimAC_header_fields)):
        if UniprotPrimAC_header_fields[i] == 'Primary_AC':
            UniProt_PrimAC_index = i
        elif UniprotPrimAC_header_fields[i] == 'TaxID':
            TaxID_index = i

    if not UniProt_PrimAC_index >= 0:
        sys.exit("Error: Missing required column title 'Primary_AC' in the file: %s \n" % inPrimAC)
    elif not TaxID_index >= 0:
        sys.exit("Error: Missing required column title 'TaxID' in the file: %s \n" % inPrimAC)
    # else grabbed the required column indices -> PROCEED

    # Data lines
    for line in UniprotPrimAC_File:
        Primary_AC_fields = line.split('\t')

        # Key -> UniProt primary accession
        # value -> Taxonomy ID
        UniprotPrimAC_dict[Primary_AC_fields[UniProt_PrimAC_index]] = Primary_AC_fields[TaxID_index]

    # Closing the file
    UniprotPrimAC_File.close()

    return UniprotPrimAC_dict

###########################################################

# Parses the tab-seperated UniProt Secondary Accession
# file produced by Uniprot_parser.py
#
# Required columns are: 'Secondary_ACs' and 'Primary_AC' (can be in any order,
# but they MUST exist)
#
# Returns a dictionary:
# - Key -> UniProt Secondary Accession
# - Value -> Corresponding UniProt Primary Accession
def SecAC(inSecAC):

    # Dictionary to store UniProt Secondary Accession and
    # UniProt Primary Accession
    UniprotSecAC_dict = {}

    UniprotSecAC_File = open(inSecAC)

    logging.info("Processing data from UniProt Secondary Accession File: %s" % inSecAC)

    # Grabbing the header line
    UniprotSecAC_header = UniprotSecAC_File.readline()

    UniprotSecAC_header = UniprotSecAC_header.rstrip('\n')

    UniprotSecAC_header_fields = UniprotSecAC_header.split('\t')

    # Check the column header and grab indexes of our columns of interest
    (UniprotSecAC_index, UniprotPrimAC_index) = (-1, -1)

    for i in range(len(UniprotSecAC_header_fields)):
        if UniprotSecAC_header_fields[i] == 'Secondary_ACs':
            UniprotSecAC_index = i
        elif UniprotSecAC_header_fields[i] == 'Primary_AC':
            UniprotPrimAC_index = i

    if not UniprotSecAC_index >= 0:
        sys.exit("Error: Missing required column title 'Secondary_ACs' in the file: %s \n" % inSecAC)
    elif not UniprotPrimAC_index >= 0:
        sys.exit("Error: Missing required column title 'Primary_AC' in the file: %s \n" % inSecAC)
    # else grabbed the required column indices -> PROCEED

    # Data lines
    for line in UniprotSecAC_File:
        line = line.rstrip("\n") # removing carriage returns
        Secondary_AC_fields = line.split('\t')

        # Key -> Uniprot Secondary_AC
        # Value -> UniProt Primary_AC
        (SecAC,PrimAC) = (Secondary_AC_fields[UniprotSecAC_index], Secondary_AC_fields[UniprotPrimAC_index])

        # If the Secondary_AC is associated with multiple Primary_AC, it is considered a bad Secondary_AC
        # No Uniprot Accession is of the type "-1"
        # So, we assign "-1" as the value to this bad Secondary_AC to avoid using it later
        if UniprotSecAC_dict.get(SecAC, False):
            if UniprotSecAC_dict[SecAC] != "-1":
                UniprotSecAC_dict[SecAC]  = "-1"
            # else: Secondary_AC is already bad => NOOP
        else:
            UniprotSecAC_dict[SecAC] = PrimAC

    # Closing the file
    UniprotSecAC_File.close()

    return  UniprotSecAC_dict

###########################################################

# Parses the tab-seperated GeneID file produced by Uniprot_parser.py
#
# Required columns are: 'GeneID' and 'Primary_AC' (can be in any order,
# but they MUST exist)
#
# Returns a dictionary:
# - Key -> Gene ID
# - Value -> Corresponding UniProt Primary Accession
def GeneID(inGeneID):

    # Dictionary to store GeneID and
    # UniProt Primary Accession
    GeneID_dict = {}

    GeneID_File = open(inGeneID)

    logging.info("Processing data from GeneID File: %s" % inGeneID)

    # Grabbing the header line
    GeneID_File_header = GeneID_File.readline()

    GeneID_File_header = GeneID_File_header.rstrip('\n')

    GeneID_File_header_fields = GeneID_File_header.split('\t')

    # Check the column header and grab indexes of our columns of interest
    (GeneID_index, UniprotPrimAC_index) = (-1, -1)

    for i in range(len(GeneID_File_header_fields)):
        if GeneID_File_header_fields[i] == 'GeneID':
            GeneID_index = i
        elif GeneID_File_header_fields[i] == 'Primary_AC':
            UniprotPrimAC_index = i

    if not GeneID_index >= 0:
        sys.exit("Error: Missing required column title 'GeneID' in the file: %s \n" % inGeneID)
    elif not UniprotPrimAC_index >= 0:
        sys.exit("Error: Missing required column title 'Primary_AC' in the file: %s \n" % inGeneID)
    # else grabbed the required column indices -> PROCEED

    # Data lines
    for line in GeneID_File:
        # removing carriage returns
        line = line.rstrip("\n")
        GeneID_fields = line.split('\t')

        # Key -> GeneID
        # Value -> Uniprot Primary Accession
        (GeneID,PrimAC) = (GeneID_fields[GeneID_index], GeneID_fields[UniprotPrimAC_index])

        # If the GeneID is associated with multiple Primary_AC, it is considered a bad GeneID
        # No GeneID is of the type "-1"
        # So, we assign "-1" as the value to this bad GeneID to avoid using it later
        if GeneID_dict.get(GeneID, False):
            if GeneID_dict[GeneID] != "-1":
                GeneID_dict[GeneID] = "-1"
            # else: GeneID is already bad => NOOP
        else:
            GeneID_dict[GeneID] = PrimAC

    # Closing the file
    GeneID_File.close()

    return GeneID_dict

###########################################################

# Parses a miTAB 2.5 or 2.7 file
# Maps the data to UniProt using the above dictionaries
#
# For each interaction, grabs the UniProt Primary Accession of the interacting Proteins
# Also grabs Interaction Detection Method, PMID and Interaction type
# Ignores interactions if:
# - UniProt Primary Accessions of either of the proteins cannot be found OR
# - if PMID is missing OR
# - if TaxID of both interacting proteins is not '9606' i.e. human
#
# Prints to STDOUT in .tsv format by sorting the Uniprot PrimAC
# The output consists of 5 columns:
# - Protein A UniProt Primary Accession
# - Protein B UniProt Primary Accession
# - Interaction Detection Method
# - Pubmed Identifier
# - Interaction Type
def interaction_parser(args):

    # Calling functions
    UniprotPrimAC_dict = PrimAC(args.inPrimAC)
    UniprotSecAC_dict = SecAC(args.inSecAC)
    GeneID_dict = GeneID(args.inGeneID)

    # Debug counters
    # To check how many times the UniProt Primary ACs are found in the respective dictionaries
    found_inPrimACFile = 0
    found_inSecACFile = 0
    found_inGeneIDFile = 0

    # When PrimAC is identified using first 2 columns of the Interaction file
    PrimAC_inMainCols = 0

    # When PrimAC is identified using AltID columns of the Interaction file
    PrimAC_inAltCols = 0

    # When PrimAC is identified using the GeneIDs
    PrimAC_foundwithGeneID = 0

    # When PrimAC is identified using the Secondary_ACs
    PrimAC_foundwithSecAC = 0

    # Keeping count of the lines where UniProt PrimAC of proteins not found
    notfound_Protein_A_PrimAC = 0
    notfound_Protein_B_PrimAC = 0

    # Keeping count of PMID not found
    notfound_PMID = 0

    # Compiling all the regular expressions###

    # uniprot ids for protein
    re_uniprot = re.compile('^uniprot(kb|/swiss-prot):([A-Z0-9-_]+)$')
    re_uniprot_missed = re.compile('^uniprot')

    # if uniprot AC not found, using GeneID to get corresponding Primary_AC
    re_GeneID = re.compile('^entrez gene/locuslink:(\d+)$')
    re_GeneID_missed = re.compile('^entrez gene/locuslink:(\d+)$')

    # PSI-MI term parser for Interaction Detection Method and Interaction type
    re_psimi = re.compile('^psi-mi:"(MI:\d+)"')
    re_psimi_missed = re.compile('^psi-mi:')

    # Pubmed Identifiers
    re_PMID = re.compile('^pubmed:(\d+)$')

    # some pubmed identifiers are unassigned in Intact (pubmed:unassigned)
    re_PMID_unassigned = re.compile('^pubmed:unassigned')
    re_PMID_missed = re.compile('^pubmed:')

    # User input Protein-protein Interaction file
    interaction_file = open(args.inInteraction)

    logging.info("Processing data from Protein-Protein Interaction File: %s" % args.inInteraction)

    # Skip header
    interaction_file.readline()

    logging.info("Preparing Output...")

    # Data lines
    for line in interaction_file:
        line = line.rstrip('\n')
        line_fields = line.split('\t')

        # Initializing variables/accumulators
        Prots = ['','']
        IntDetectMethod = ''
        PMID = ''
        Interaction_type = ''

        # line_fields[15] -> Expansion method(s) in mitab 2.7
        # For true binary interactions, this column will be '-'
        # thus, eliminates spoke expansion
        if (len(line_fields)>15) and (line_fields[15] != '-'):
            continue
        # else if mitab 2.7 no spoke expansion
        for protindex in [0,1]:
            if (re_uniprot.match(line_fields[protindex])):
                ID = re_uniprot.match(line_fields[protindex]).group(2)
                # Check if it exists in the dictionary
                if UniprotPrimAC_dict.get(ID, False):
                    Prots[protindex] = ID
                    found_inPrimACFile += 1
                    PrimAC_inMainCols += 1
                    continue
            elif (re_uniprot_missed.match(line_fields[protindex])):
                sys.exit("Error: ID is a uniprot Accession but failed to grab it for the line:\n" + line)
            elif (re_GeneID.match(line_fields[protindex])):
                ID = re_GeneID.match(line_fields[protindex]).group(1)
                # Check if it exists in the dictionary and isn't bad ie "-1"
                if GeneID_dict.get(ID, "-1") != "-1":
                    Prots[protindex] = GeneID_dict[ID]
                    found_inGeneIDFile += 1
                    PrimAC_inMainCols += 1
                    PrimAC_foundwithGeneID += 1
                    continue
            elif (re_GeneID_missed.match(line_fields[protindex])):
                sys.exit("Error: ID is a GeneID but failed to grab it for the line:\n" + line)

            # Uniprot AC not found/not primary_AC and GeneID not found,
            # then, look in Alternate ID columns of the interaction_file
            # Stop when the first AltID is recognized (Careful, if multiple altIDs would be found)
            altIDs = line_fields[2+protindex].split('|')
            for altID in altIDs:
                if (re_uniprot.match(altID)):
                    ID = re_uniprot.match(altID).group(2)
                    # Check if it exists in the dictionary
                    if UniprotPrimAC_dict.get(ID, False):
                        Prots[protindex] = ID
                        found_inPrimACFile += 1
                        PrimAC_inAltCols += 1
                        # we want "next protindex" but python doesn't have this
                        # So we break and exit the loop
                        break
                    # ElseIf the accession is found in the Secondary_AC_dict
                elif UniprotSecAC_dict.get(ID, "-1") != "-1":
                        # Use the corresponding Primary AC
                        Prots[protindex] = UniprotSecAC_dict[ID]
                        found_inSecACFile += 1
                        PrimAC_inAltCols += 1
                        PrimAC_foundwithSecAC += 1
                        break
                elif (re_uniprot_missed.match(altID)):
                    sys.exit("Error: AltID "+altID+" is Uniprot Accession but failed to grab it for line:\n" + line)
                elif (re_GeneID.match(altID)):
                    ID = re_GeneID.match(altID).group(1)
                    # Check if it exists in the dictionary and isn't bad ie "-1"
                    if GeneID_dict.get(ID, "-1") != "-1":
                        Prots[protindex] = GeneID_dict[ID]
                        found_inGeneIDFile += 1
                        PrimAC_inAltCols += 1
                        PrimAC_foundwithGeneID += 1
                        break
                elif (re_GeneID_missed.match(altID)):
                    sys.exit("Error: AltID "+altID+" is a GeneID but failed to grab it for line:\n" + line)
                # else: altID not recognized, look at next altID ie NOOP

        if Prots[0] == '':
            notfound_Protein_A_PrimAC += 1 # Keep count of the missing UniProt PrimAC of protein A
        if Prots[1] == '':
            notfound_Protein_B_PrimAC += 1 # Keep count of the missing UniProt PrimAC of protein B

        # if either Uniprot PrimAC is not found
        if Prots[0] == '' or Prots[1] == '':
            continue # to next line

        # if we get here both partners were found, grab the remaining data

        # Grab Interaction Detection Method
        if re_psimi.match(line_fields[6]):
            IntDetectMethod = re_psimi.match(line_fields[6]).group(1)
        elif re_psimi_missed.match(line_fields[6]):
            sys.exit("Error: Failed to grab the Interaction Detection Method for line:\n" + line)

        # Grab Pubmed Identifier
        # Some Publication Identifier fields include additional fields such as MINT
        # Ex: pubmed:10542231|mint:MINT-5211933, so we split at '|'
        PMID_fields = line_fields[8].split('|')
        for PMID_entry in PMID_fields:
            if (re_PMID.match(PMID_entry)):
                PMID = re_PMID.match(PMID_entry).group(1)
            elif (re_PMID_unassigned.match(PMID_entry)):
                continue
            elif (re_PMID_missed.match(PMID_entry)):
                sys.exit("Error: Failed to grab the Pubmed Id for line:\n" + line)

        # Grab Interaction type
        if (re_psimi.match(line_fields[11])):
            Interaction_type = re_psimi.match(line_fields[11]).group(1)
        elif (re_psimi_missed.match(line_fields[11])):
            sys.exit("Error: Failed to grab the Interaction_type for line:\n" + line)

        # If no PMID, skip line
        if (PMID == ''):
            notfound_PMID += 1
            continue

        # We want only Human-Human Interaction experiments
        if (UniprotPrimAC_dict[Prots[0]] != '9606') or (UniprotPrimAC_dict[Prots[1]] != '9606'):
            continue

        # Here we grabbed all the necessary data, print to STDOUT and move on to the next line
        if Prots[0] < Prots[1]: # Sorting the PrimAC
            interaction_out_line = [Prots[0], Prots[1], IntDetectMethod, PMID, Interaction_type]
        else:
            interaction_out_line = [Prots[1], Prots[0], IntDetectMethod, PMID, Interaction_type]

        print("\t".join(interaction_out_line))

    # Debug counters
    logging.debug("No. of times Uniprot Primary Accession not found for Protein A: %d" % notfound_Protein_A_PrimAC)
    logging.debug("No. of times Uniprot Primary Accession not found for Protein B: %d" % notfound_Protein_B_PrimAC)
    logging.debug("No. of times Pubmed ID not found for the Interaction: %d" % notfound_PMID)
    logging.debug("No. of times Uniprot Primary Accession identified using the Uniprot Primary Accession File: %d" % found_inPrimACFile)
    logging.debug("No. of times Uniprot Primary Accession identified using the Uniprot Secondary Accession File: %d" % found_inSecACFile)
    logging.debug("No. of times Uniprot Primary Accession identified using the GeneID File: %d" % found_inGeneIDFile)
    logging.debug("No. of times Uniprot Primary Accession identified using the first 2 columns of the Interaction file: %d" % PrimAC_inMainCols)
    logging.debug("No. of times Uniprot Primary Accession identified using the AltID columns of the Interaction file: %d" % PrimAC_inAltCols)
    logging.debug("No. of times Uniprot Primary Accession identified using GeneIDs: %d" % PrimAC_foundwithGeneID)
    logging.debug("No. of times Uniprot Primary Accession identified using Secondary_AC: %d" % PrimAC_foundwithSecAC)

    # Closing the file
    interaction_file.close()

    logging.info("Done 🎉")

    return

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
-------------------------------------------------------------------------------------------------------
Program: Parses a miTAB 2.5 or 2.7 file, maps to the uniprot file and prints to STDOUT in tsv format
-------------------------------------------------------------------------------------------------------
The output (Human-Human Protein Interaction Experiments) consists of five columns in .tsv format:
 -> UniProt Primary Accession of Protein A
 -> UniProt Primary Accession of Protein B
 -> Interaction Detection Method
 -> Pubmed Identifier
 -> Interaction type
 -----------------------------------------------------------------------------------------------------

 Arguments [defaults] -> Can be abbreviated to shortest unambiguous prefixes
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inInteraction', metavar = "Input File", dest = "inInteraction", help = 'Input File Name (Protein-protein Interaction file)', required = True)
    required.add_argument('--inPrimAC', metavar = "Input File", dest = "inPrimAC", help = 'Uniprot Primary Accession File generated by the uniprot parser', required = True)
    required.add_argument('--inSecAC', metavar = "Input File", dest = "inSecAC", help = 'Uniprot Secondary Accession File generated by the uniprot parser', required = True)
    required.add_argument('--inGeneID', metavar = "Input File", dest = "inGeneID", help = 'GeneID File generated by the uniprot parser', required = True)

    args = file_parser.parse_args()
    interaction_parser(args)

if __name__ == "__main__":
    # Logging to Standard Error
    Log_Format = "%(levelname)s - %(asctime)s - %(message)s \n"
    logging.basicConfig(stream = sys.stderr, format  = Log_Format, level = logging.DEBUG)
    main()
