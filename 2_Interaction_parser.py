#!/usr/bin/python

# manojmw
# 14 Jan, 2022

import re
import sys, argparse
import logging

###########################################################

# Parses the tab-seperated UniProt Output file
# produced by 1_Uniprot_parser.py
#
# Required columns are: 
# - 'Primary_AC' 
# - 'TaxID' 
# - 'Secondary_ACs'
# - 'GeneIDs'
# - 'GeneNames'
# (can be in any order,
# but they MUST exist)
#
# Returns 4 dictionaries:
# UniprotPrimAC_dict dictionary contains:
# - Key -> UniProt Primary Accession
# - Value -> TaxID
#
# UniprotSecAC_dict dictionary contains:
# - Key -> UniProt Secondary Accession
# - Value -> Corresponding UniProt Primary Accession
# 
# GeneID_dict dictionary contains:
# - Key -> GeneID
# - Value -> Corresponding UniProt Primary Accession
#
# GeneName_dict dictionary contains:
# - Key -> GeneName
# - Value -> Corresponding UniProt Primary Accession
def Build_UniProtDicts(inUniProt):

    logging.info("Starting to run...")

    Uniprot_File = open(inUniProt)

    logging.info("Processing data from UniProt File: %s" % inUniProt)

    # Grabbing the header line
    Uniprot_header = Uniprot_File.readline()

    Uniprot_header = Uniprot_header.rstrip('\n')

    Uniprot_header_fields = Uniprot_header.split('\t')

    # Dictionary to store UniProt Primary Accession
    # and the Taxonomy Identifier
    UniprotPrimAC_dict = {}

    # Dictionary to store UniProt Secondary Accession and
    # UniProt Primary Accession
    UniprotSecAC_dict = {}

    # Dictionary to store GeneID and
    # UniProt Primary Accession
    GeneID_dict = {}

    # Dictionary to store GeneName and
    # UniProt Primary Accession
    GeneName_dict = {}

    (UniProtPrimAC_index, TaxID_index, UniprotSecAC_index, GeneID_index, GeneName_index) = (-1, -1, -1, -1, -1)

    # Check the column header and grab indexes of our columns of interest
    for i in range(len(Uniprot_header_fields)):
        if Uniprot_header_fields[i] == 'Primary_AC':
            UniProtPrimAC_index = i
        elif Uniprot_header_fields[i] == 'TaxID':
            TaxID_index = i
        elif Uniprot_header_fields[i] == 'Secondary_ACs':
            UniprotSecAC_index = i
        elif Uniprot_header_fields[i] == 'GeneIDs':
            GeneID_index = i
        elif Uniprot_header_fields[i] == 'GeneNames':
            GeneName_index = i    

    # Sanity check
    if not UniProtPrimAC_index >= 0:
        sys.exit("Error: Missing required column title 'Primary_AC' in the file: %s \n" % inUniProt)
    elif not TaxID_index >= 0:
        sys.exit("Error: Missing required column title 'TaxID' in the file: %s \n" % inUniProt)
    elif not UniprotSecAC_index >= 0:
        sys.exit("Error: Missing required column title 'Secondary_ACs' in the file: %s \n" % inUniProt)
    elif not GeneID_index >= 0:
        sys.exit("Error: Missing required column title 'GeneIDs' in the file: %s \n" % inUniProt)
    elif not GeneName_index >= 0:
        sys.exit("Error: Missing required column title 'GeneNames' in the file: %s \n" % inUniProt)
    # else grabbed the required column indices -> PROCEED

    # Data lines
    for line in Uniprot_File:
        UniProt_line = line.rstrip("\n") # removing carriage returns
        UniProt_fields = UniProt_line.split('\t')

        # Populating UniprotPrimAC_dict
        # Key -> UniProt primary accession
        # value -> Taxonomy ID
        UniprotPrimAC_dict[UniProt_fields[UniProtPrimAC_index]] = UniProt_fields[TaxID_index]

        # Populating UniprotSecAC_dict
        # Key -> Uniprot Secondary_AC
        # Value -> UniProt Primary_AC
        # Secondary_ACs field in the inFile can contain 
        # a single UniProt Secondary Accession or a comma-seperated
        # list of UniProt Secondary Accessions
        try:
            UniProt_SecACs = UniProt_fields[UniprotSecAC_index].split(',')
        except:
            UniProt_SecACs = [UniProt_fields[UniprotSecAC_index]]

        for UniProt_SecAC in UniProt_SecACs:
            (SecAC,PrimAC) = (UniProt_SecAC, UniProt_fields[UniProtPrimAC_index])

            # If the UniProt Secondary_AC is associated with multiple Primary_AC, 
            # it is considered a bad Secondary_AC 
            # No Uniprot Accession is of the type "-1"
            # So, we assign "-1" as the value to this bad Secondary_AC to avoid using it later
            if UniprotSecAC_dict.get(SecAC, False):
                if UniprotSecAC_dict[SecAC] != "-1":
                    UniprotSecAC_dict[SecAC]  = "-1"
                # else: Secondary_AC is already bad => NOOP
            else:
                UniprotSecAC_dict[SecAC] = PrimAC

        # Populating GeneID_dict
        # Key -> GeneID
        # Value -> Uniprot Primary Accession
        # GeneIDs field in the inFile can contain 
        # a single GeneID or a comma-seperated
        # list of GeneIDs        
        try:
            UniProt_GeneIDs = UniProt_fields[GeneID_index].split(',')
        except:
            UniProt_GeneIDs = [UniProt_fields[GeneID_index]]

        for UniProt_GeneID in UniProt_GeneIDs:
            (GeneID,PrimAC) = (UniProt_GeneID, UniProt_fields[UniProtPrimAC_index])

            # If the GeneID is associated with multiple Primary_AC, 
            # it is considered a bad GeneID
            # No GeneID is of the type "-1"
            # So, we assign "-1" as the value to this bad GeneID to avoid using it later
            if GeneID_dict.get(GeneID, False):
                if GeneID_dict[GeneID] != "-1":
                    GeneID_dict[GeneID] = "-1"
                # else: GeneID is already bad => NOOP
            else:
                GeneID_dict[GeneID] = PrimAC

        # Populating GeneName_dict
        # Key -> GeneName
        # Value -> Uniprot Primary Accession
        # GeneIDs field in the inFile can contain 
        # a single GeneID or a comma-seperated
        # list of GeneIDs        
        try:
            UniProt_GeneNames = UniProt_fields[GeneName_index].split(',')
        except:
            UniProt_GeneNames = [UniProt_fields[GeneName_index]]

        for UniProt_GeneName in UniProt_GeneNames:
            (GeneName,PrimAC) = (UniProt_GeneName, UniProt_fields[UniProtPrimAC_index])

            # If the GeneName is associated with multiple Primary_AC, 
            # it is considered a bad GeneName
            # No GeneName is of the type "-1"
            # So, we assign "-1" as the value to this bad GeneName to avoid using it later
            if GeneName_dict.get(GeneName, False):
                if GeneName_dict[GeneName] != "-1":
                    GeneID_dict[GeneName] = "-1"
                # else: GeneID is already bad => NOOP
            else:
                GeneName_dict[GeneName] = PrimAC        

    # Closing the file
    Uniprot_File.close()

    return UniprotPrimAC_dict, UniprotSecAC_dict, GeneID_dict, GeneName_dict


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
    (UniprotPrimAC_dict, UniprotSecAC_dict, GeneID_dict, GeneName_dict) = Build_UniProtDicts(args.inUniProt)

    # When PrimAC is identified using first 2 columns of the Interaction file
    PrimAC_inMainCols = 0

    # When PrimAC is identified using AltID columns of the Interaction file
    PrimAC_inAltCols = 0

    # When PrimAC is identified using AltID columns of the Interaction file
    PrimAC_inAliasCols = 0

    # When PrimAC is identified using the GeneIDs
    PrimAC_foundwithGeneID = 0

    # When PrimAC is identified using the Secondary_ACs
    PrimAC_foundwithSecAC = 0

    # When PrimAC is identified using the Gene Name
    PrimAC_foundwithGeneName = 0

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

    # Getting UniProt Primary AC using Gene Name if not found using above
    re_GeneName = re.compile('^(entrez gene.locuslink:|uniprotkb:)([\w\s\_\\\/\:\.\-]+$|[\w\s\_\\\/\:\.\-]+)') 

    # PSI-MI term parser for Interaction Detection Method and Interaction type
    re_psimi = re.compile('^psi-mi:"(MI:\d+)"')
    re_psimi_missed = re.compile('^psi-mi:')

    # Pubmed Identifiers
    re_PMID = re.compile('^pubmed:(\d+)$')

    # some pubmed identifiers are unassigned in Intact (pubmed:unassigned)
    re_PMID_unassigned = re.compile('^pubmed:unassigned')
    re_PMID_missed = re.compile('^pubmed:')

    # User input Protein-protein Interaction file
    interaction_file = open(args.inInteraction, encoding="utf-8")

    logging.info("Processing data from Protein-Protein Interaction File: %s" % args.inInteraction)

    # Skip header
    interaction_file.readline()

    # Data lines
    for line in interaction_file:
        line = line.rstrip('\n')
        line_fields = line.split('\t')

        # Initializing variables/accumulators
        Prots = ['','']
        IntDetectMethod = ''
        PMID = ''
        Interaction_type = ''

        # Using 'Complex expansion' column (for miTAB 2.7) to eliminate expansion data
        # is not reliable and can sometimes result in the loss of true binary interactions
        # as observed in miTAB 2.7 file retrieved from IntAct (at the time of
        # writing this script. see below example)
        # Example for the interactions between proteins (from IntAct): 
        # P81274 & Q14980 
        # Interaction Detection Method: 2 hybrid 
        # Pubmed:15537540 
        # Interaction AC: EBI-624047 
        # Expansion: spoke expansion
        # Interaction Detection Method is 2-Hybrid and this is marked as
        # expansion which is most probably incorrect
        #
        # In some cases, Affinity Chromatography Technology (ACT) is not marked as
        # expansion, although most of the PPIs determined by ACT 
        # correspond to expansion
        # Example (from IntAct):
        # P19784 & P03211 
        # Interaction Detection Method: Affinity Chromatography Technology
        # Pubmed: 12783858 
        # Interaction AC: EBI-8656048
        # Expansion: Not spoke expansion
        #
        # So through private communication with a PPI database,
        # we have set a few criteria to eliminate the expansion data while
        # buliding the Interactome in the Build_Interactome.py script
        for protindex in [0,1]:
            if (re_uniprot.match(line_fields[protindex])):
                ID = re_uniprot.match(line_fields[protindex]).group(2)
                # Check if it exists in the dictionary
                if UniprotPrimAC_dict.get(ID, False):
                    Prots[protindex] = ID
                    PrimAC_inMainCols += 1
                    continue
            elif (re_uniprot_missed.match(line_fields[protindex])):
                sys.exit("Error: ID is a uniprot Accession but failed to grab it for the line:\n" + line)
            elif (re_GeneID.match(line_fields[protindex])):
                ID = re_GeneID.match(line_fields[protindex]).group(1)
                # Check if it exists in the dictionary and isn't bad ie "-1"
                if GeneID_dict.get(ID, "-1") != "-1":
                    Prots[protindex] = GeneID_dict[ID]
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
                        PrimAC_inAltCols += 1
                        # we want "next protindex" but python doesn't have this
                        # So we break and exit the loop
                        break
                        # ElseIf the accession is found in the Secondary_AC_dict
                    elif UniprotSecAC_dict.get(ID, "-1") != "-1":
                        # Use the corresponding Primary AC
                        Prots[protindex] = UniprotSecAC_dict[ID]
                        PrimAC_inAltCols += 1
                        PrimAC_foundwithSecAC += 1
                        break
                elif (re_uniprot_missed.match(altID)):
                    sys.exit("Error: AltID "+altID+" is Uniprot Accession but failed to grab it for the line:\n" + line)
                elif (re_GeneID.match(altID)):
                    GeneID = re_GeneID.match(altID).group(1)
                    # Check if it exists in the GeneID dictionary and isn't bad ie "-1"
                    if GeneID_dict.get(GeneID, "-1") != "-1":
                        Prots[protindex] = GeneID_dict[GeneID]
                        PrimAC_inAltCols += 1
                        PrimAC_foundwithGeneID += 1
                        break
                elif (re_GeneID_missed.match(altID)):
                    sys.exit("Error: AltID "+altID+" is a GeneID but failed to grab it for the line:\n" + line)
                    # If UniProt PrimAC not found using above,
                    # then find UniProt PrimAC using the Gene name    
                    # In case of miTAB 2.5 format, Gene name can be
                    # present in AltIDs column as well (as seen in BioGRID)
                    # Where as in miTAB 2.7 (as seen in IntAct), both 
                    # Gene name and Gene name Synonym is present in Alias column
                    # So we check in AltID columns
                elif(re_GeneName.match(altID)):
                    GN_altID = re_GeneName.match(altID).group(2)
                    # Check if it exists in the Gene Name dictionary and isn't bad ie "-1"
                    if GeneName_dict.get(GN_altID, "-1") != "-1":
                        Prots[protindex] = GeneName_dict[GN_altID]
                        PrimAC_inAltCols += 1
                        PrimAC_foundwithGeneName += 1
                        break 
                # else: altID not recognized, look at next altID ie NOOP

            # If UniProt Primary AC not found in the AltID column either, 
            # then look in the Alias(es) column
            aliasIDs = line_fields[4+protindex].split('|')
            for aliasID in aliasIDs:  
                if (re_GeneName.match(aliasID)):
                    GN_aliasID = re_GeneName.match(aliasID).group(2)
                    # Check if it exists in the Gene Name dictionary and isn't bad ie "-1"
                    if GeneName_dict.get(GN_aliasID, "-1") != "-1":
                        Prots[protindex] = GeneName_dict[GN_aliasID]
                        PrimAC_inAliasCols += 1
                        PrimAC_foundwithGeneName += 1
                        break 

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
                sys.exit("Error: Failed to grab the Pubmed Id for the line:\n" + line)

        # Grab Interaction type
        if (re_psimi.match(line_fields[11])):
            Interaction_type = re_psimi.match(line_fields[11]).group(1)
        elif (re_psimi_missed.match(line_fields[11])):
            sys.exit("Error: Failed to grab the Interaction_type for the line:\n" + line)

        # If no PMID, skip line
        if (PMID == ''):
            notfound_PMID += 1
            continue

        # Sanity check: only Human-Human Interaction experiments
        if (UniprotPrimAC_dict[Prots[0]] != '9606') or (UniprotPrimAC_dict[Prots[1]] != '9606'):
            continue

        # Here we grabbed all the necessary data, print to STDOUT and move on to the next line
        if Prots[0] < Prots[1]: # Sorting the PrimAC
            interaction_out_line = [Prots[0], Prots[1], IntDetectMethod, PMID, Interaction_type]
        else:
            interaction_out_line = [Prots[1], Prots[0], IntDetectMethod, PMID, Interaction_type]

        print("\t".join(interaction_out_line))

    # Debug counters
    # logging.debug("No. of times Uniprot Primary Accession not found for Protein A: %d" % notfound_Protein_A_PrimAC)
    # logging.debug("No. of times Uniprot Primary Accession not found for Protein B: %d" % notfound_Protein_B_PrimAC)
    # logging.debug("No. of times Pubmed ID not found for the Interaction: %d" % notfound_PMID)
    # logging.debug("No. of times Uniprot Primary Accession identified using the first 2 columns of the Interaction file: %d" % PrimAC_inMainCols)
    # logging.debug("No. of times Uniprot Primary Accession identified using the AltID columns of the Interaction file: %d" % PrimAC_inAltCols)
    # logging.debug("No. of times Uniprot Primary Accession identified using the Alias(es) columns of the Interaction file: %d" % PrimAC_inAliasCols)
    # logging.debug("No. of times Uniprot Primary Accession identified using GeneIDs: %d" % PrimAC_foundwithGeneID)
    # logging.debug("No. of times Uniprot Primary Accession identified using Secondary_AC: %d" % PrimAC_foundwithSecAC)
    # logging.debug("No. of times Uniprot Primary Accession identified using Gene Name: %d" % PrimAC_foundwithGeneName)

    # Closing the file
    interaction_file.close()

    logging.info("All Done, completed successfully!")

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

    required.add_argument('--inInteraction', metavar = "Input File", dest = "inInteraction", help = 'Input File Name (Protein-protein Interaction file: miTAB 2.5 OR 2.7)', required = True)
    required.add_argument('--inUniProt', metavar = "Input File", dest = "inUniProt", help = 'Uniprot Output File generated by the Uniprot_parser.py', required = True)

    args = file_parser.parse_args()
    interaction_parser(args)

if __name__ == "__main__":
    # Logging to Standard Error
    Log_Format = "%(levelname)s - %(asctime)s - %(message)s \n"
    logging.basicConfig(stream = sys.stderr, format  = Log_Format, level = logging.DEBUG)
    main()
