#!/usr/bin/python

import re
import argparse
import sys

###Creating dictionaries from uniprot output files###

###UniProt Primary Accession dictionary
def PrimAC(args):
    Primary_AC_dict = {} ##Initializing an empty dictionary
    Primary_AC_file = open(args.inPrimAC)

    Primary_AC_file.readline() ###Skip header

    for line in Primary_AC_file:
        Primary_AC_fields = line.split('\t')
        (PrimAC,TaxID) = (Primary_AC_fields[0], Primary_AC_fields[1]) ##Key -> primary accession
        Primary_AC_dict[PrimAC] = TaxID

    return Primary_AC_dict


###UniProt Secondary Accession dictionary
def SecAC(args):
    Secondary_AC_dict = {} ##Initializing an empty dictionary
    Secondary_AC_file = open(args.inSecAC)

    Secondary_AC_file.readline() ###Skip header

    for line in Secondary_AC_file:
        line = line.rstrip("\n") ##removing carriage returns
        Secondary_AC_fields = line.split('\t')
        (SecAC,PrimAC) = (Secondary_AC_fields[1], Secondary_AC_fields[0]) ##Key -> secondary accession
        ###If the Secondary_AC is associated with multiple Primary_AC, it is considered a bad Secondary_AC
        ###No Uniprot Accession is of the type "-1"
        ###So, we assign "-1" as the value to this bad Secondary_AC to prevent it from being added to the Secondary_AC_dict
        if Secondary_AC_dict.get(SecAC, False):
            if Secondary_AC_dict[SecAC] != "-1":
                Secondary_AC_dict[SecAC]  = "-1"
            #else: Secondary_AC is already bad => NOOP
        else:
            Secondary_AC_dict[SecAC] = PrimAC
    return  Secondary_AC_dict

###GeneID dictionary
def GeneID(args):
    GeneID_dict = {} ##Initializing an empty dictionary
    GeneID_file = open(args.inGeneID)

    GeneID_file.readline() ###Skip header

    for line in GeneID_file:
        line = line.rstrip("\n") ##removing carriage returns
        GeneID_fields = line.split('\t')
        (GeneID,PrimAC) = (GeneID_fields[1], GeneID_fields[0]) ##Key -> GeneID
                ###If the GeneID is associated with multiple Primary_AC, it is considered a bad GeneID
                ###No GeneID is of the type "-1"
                ###So, we assign "-1" as the value to this bad GeneID to prevent it from being added to the GeneID_dict
        if GeneID_dict.get(GeneID, False):
            if GeneID_dict[GeneID] != "-1":
                GeneID_dict[GeneID] = "-1"
            #else: GeneID is already bad => NOOP
        else:
            GeneID_dict[GeneID] = PrimAC
    return GeneID_dict

###Protein-Protein Interaction Parser
def interaction_parser(args):
    try:
        with open(args.output_interaction, 'w') as int_out_file, open(args.output_debugcount, 'w') as debug_out_file:
            header_interaction = ['Protein_A_UniprotPrimAC', 'Protein_B_UniprotPrimAC', 'Interaction_Detect_Method', 'PMID', 'Interaction_type']
            print("\t".join(header_interaction), file = int_out_file)

            ###Calling dictionary functions
            Primary_AC_dict = PrimAC(args)
            Secondary_AC_dict = SecAC(args)
            GeneID_dict = GeneID(args)

            ###Debug counters
            ###To check how many times the ACs are found in the respective dictionaries
            found_inPrimACFile = 0
            found_inSecACFile = 0
            found_inGeneIDFile = 0

            ###User input -> Protein-protein Interaction file
            interaction_file = open(args.inInteraction)

            ###Skip header
            interaction_file.readline()

            ###Compiling all the regular expressions###

            ###uniprot ids for protein
            re_uniprot = re.compile('^uniprot(kb|/swiss-prot):([A-Z0-9-_]+)$')
            re_uniprot_missed = re.compile('^uniprot')
            ###if uniprot AC not found, using GeneID to get corresponding Primary_AC
            re_GeneID = re.compile('^entrez gene/locuslink:(\d+)$')
            re_GeneID_missed = re.compile('^entrez gene/locuslink:(\d+)$')
            ###PSI-MI term parser for Interaction Detection Method and Interaction type
            re_psimi = re.compile('^psi-mi:"(MI:\d+)"')
            re_psimi_missed = re.compile('^psi-mi:')
            ###Pubmed Identifiers
            re_PMID = re.compile('^pubmed:(\d+)$')
            ###some pubmed identifiers are unassigned (pubmed:unassigned)
            re_PMID_missed = re.compile('^pubmed:^(unassigned)')

            ###Parsing the interaction file
            for line in interaction_file:
                line = line.rstrip('\n')
                line_fields = line.split('\t')

                ###Initializing variables/accumulators
                Prots = ['','']
                IntDetectMethod = ''
                PMID = ''
                Interaction_type = ''

                for protindex in [0,1]:
                    if (re_uniprot.match(line_fields[protindex])):
                        ID = re_uniprot.match(line_fields[protindex]).group(2)
                        ##Check if it exists in the dictionary
                        if Primary_AC_dict.get(ID, False):
                            Prots[protindex] = ID
                            found_inPrimACFile = found_inPrimACFile + 1
                            continue
                    elif (re_uniprot_missed.match(line_fields[protindex])):
                        sys.exit("ID is a uniprot Accession but failed to grab it for the line:\n" + line)
                    elif (re_GeneID.match(line_fields[protindex])):
                        ID = re_GeneID.match(line_fields[protindex]).group(1)
                        ##Check if it exists in the dictionary
                        if GeneID_dict.get(ID, False):
                            Prots[protindex] = GeneID_dict[ID]
                            found_inGeneIDFile = found_inGeneIDFile + 1
                            continue
                    elif (re_GeneID_missed.match(line_fields[protindex])):
                        sys.exit("ID is a GeneID but failed to grab it for the line:\n" + line)

                    ###Uniprot AC not found/not primary_AC and GeneID not found,
                    ###then, look in Alternate ID columns of the interaction_file
                    altIDs = line_fields[2+protindex].split('|')
                    for altID in altIDs:
                        if (re_uniprot.match(altID)):
                            ID = re_uniprot.match(altID).group(2)
                            ##Check if it exists in the dictionary
                            if Primary_AC_dict.get(ID, False):
                                Prots[protindex] = ID
                                found_inPrimACFile = found_inPrimACFile + 1
                                # we want "next protindex" but python doesn't have this
                                ##So we break and exit the loop
                                break
                            ###ElseIf the accession is found in the Secondary_AC_dict
                            elif Secondary_AC_dict.get(ID, "-1") != "-1":
                                ###Use the corresponding Primary AC
                                Prots[protindex] = Secondary_AC_dict[ID]
                                found_inSecACFile = found_inSecACFile + 1
                                break
                        elif (re_uniprot_missed.match(altID)):
                            sys.exit("altID "+altID+" is Uniprot Accession but failed to grab it for line:\n" + line)
                        elif (re_GeneID.match(altID)):
                            ID = re_GeneID.match(altID).group(1)
                            ##Check if it exists in the dictionary
                            if GeneID_dict.get(ID, False):
                                Prots[protindex] = GeneID_dict[ID]
                                found_inGeneIDFile = found_inGeneIDFile + 1
                                break
                        elif (re_GeneID_missed.match(altID)):
                            sys.exit("altID "+altID+" is a GeneID but failed to grab it for line:\n" + line)
                        #else: altID not recognized, look at next altID ie NOOP
                    #if protindex==0 and we didn't recognize the partner, no point in looking for second partner
                    if Prots[protindex] == '':
                        break
                #done looking for both partners, if we didn't find them both move to next line
                if Prots[1] == '':
                    #testing Prots[1] is sufficient, because if Prots[0] not found we broke
                    continue #to next line

                #if we get here both partners were found, grab the remaining data
                if re_psimi.match(line_fields[6]):
                    IntDetectMethod = re_psimi.match(line_fields[6]).group(1)
                elif re_psimi_missed.match(line_fields[6]):
                    sys.exit("Failed to grab the Interaction Detection Method for line:\n" + line)
                ###Some Publication Identifier fields include additional fields such as MINT
                ###Ex: pubmed:10542231|mint:MINT-5211933, so we split at '|'
                PMID_fields = line_fields[8].split('|')
                for PMIDs in PMID_fields:
                    if (re_PMID.match(PMIDs)):
                        PMID = re_PMID.match(PMIDs).group(1)
                    elif (re_PMID_missed.match(PMIDs)):
                        sys.exit("Failed to grab the Pubmed Id for line:\n" + line)
                if (re_psimi.match(line_fields[11])):
                    Interaction_type = re_psimi.match(line_fields[11]).group(1)
                elif (re_psimi_missed.match(line_fields[11])):
                    sys.exit("Failed to grab the Interaction_type for line:\n" + line)

                #Here we grabbed all the necessary data, print to the file and move on to the next line
                interaction_out_line = [Prots[0], Prots[1], IntDetectMethod, PMID, Interaction_type]
                print("\t".join(interaction_out_line) , file = int_out_file)

            ###Debug counter
            print("No. of times Uniprot Primary Accession found in the Uniprot Primary Accession File:", found_inPrimACFile, file = debug_out_file)
            print("No. of times Uniprot Primary Accession found in the Uniprot Secondary Accession File:", found_inSecACFile, file = debug_out_file)
            print("No. of times Uniprot Primary Accession found in the GeneID File:", found_inGeneIDFile, file = debug_out_file)

        int_out_file.close()

    except IOError as E:
        print("Error: Unable to open the file for writing")


####Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
----------------------------------------------------------------------------------------------
Program: Parses a MITAB 2.5 or 2.7 file, maps to the uniprot file and produces 2 output files
----------------------------------------------------------------------------------------------
Output File 1 (--outInteraction): A tab-seperated file (.tsv) with five columns
                                   -> Protein A UniProt Primary Accession
                                   -> Protein B UniProt Primary Accession
                                   -> Interaction Detection Method
                                   -> Pubmed Identifier
                                   -> Interaction type

Output File 2 (--outCount):       This file gives the count of how often the UniProt Primary
                                  Accessions are found in input mapping files i.e.
                                   -> Uniprot Primary Accession File or
                                   -> Uniprot Secondary Accession File or
                                   -> GeneID File
----------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inInteraction', metavar = "Input File", dest = "inInteraction", help = 'Input File Name (Protein-protein Interaction file)', required = True)
    required.add_argument('--inPrimaryAC', metavar = "Input File", dest = "inPrimAC", help = 'Uniprot Primary Accession File generated by the uniprot parser', required = True)
    required.add_argument('--inSecondaryAC', metavar = "Input File", dest = "inSecAC", help = 'Uniprot Secondary Accession File generated by the uniprot parser', required = True)
    required.add_argument('--inGeneID', metavar = "Input File", dest = "inGeneID", help = 'GeneID File generated by the uniprot parser', required = True)
    required.add_argument('--outInteraction', metavar = "Output File", dest = "output_interaction", help = 'Output File containing the interactions mapped to Uniprot files', required = True)
    required.add_argument('--outCount', metavar = "Output File", dest = "output_debugcount", help = 'Output File containing frequency of occurence of Uniprot Accession', required = True)
    file_parser.set_defaults(func=interaction_parser)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
