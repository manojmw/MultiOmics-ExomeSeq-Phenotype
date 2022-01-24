#!/usr/bin/python

import re
import csv
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

    ###Calling dictionary functions
    Primary_AC_dict = PrimAC(args)
    Secondary_AC_dict = SecAC(args)
    GeneID_dict = GeneID(args)

    ###User input -> Protein-protein Interaction file
    interaction_file = open(args.inInteraction)

    ###Skip header
    interaction_file.readline()

    ###Compiling all the regular expressions###

    ###uniprot ids for protein
    re_uniprot = re.compile('^uniprot(kb|/swiss-prot):([A-Z0-9]+)$')
    re_uniprot_missed = re.compile('^uniprot')
    re_GeneID = re.compile('^entrez gene/locuslink:(\d+)$')
    re_GeneID_missed = re.compile('^entrez gene/locuslink')
    re_IntDetectMethod = re.compile('^psi-mi:"(MI:\d+)$"')
    re_psimi_missed = re.compile('^psi-mi:')
    re_PMID = re.compile('^pubmed:(\d+)$')
    re_PMID_missed = re.compile('^pubmed:')
    re_IntType = re.compile('^psi-mi:"(MI:\d+)$"')

    for line in interaction_file:
        line = line.rstrip('\n')
        line_fields = line.split('\t')

        ###Initializing variables/accumulators
        Prots = ['','']
        IntDetectMethod = ''
        PMID = 0
        Interaction_type = ''

        break_flag = False

        for protindex in [0,1]:
            if (re_uniprot.match(line_fields[protindex])):
                ID = re_uniprot.match(line_fields[protindex]).group(2)
                ##Check if it exists in the dictionary
                if Primary_AC_dict.get(ID, False):
                    Prots[protindex] = ID
                    continue
            elif (re_uniprot_missed.match(line_fields[protindex])):
                sys.exit("ID is a uniprot Accession but failed to grab it for the line:\n" + line)
            elif (re_GeneID.match(line_fields[protindex])):
                ID = re_GeneID.match(line_fields[protindex]).group(1)
                ##Check if it exists in the dictionary
                if GeneID_dict.get(ID, False):
                    Prots[protindex] = GeneID_dict[ID]
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
                        if Prots[protindex] == '':
                            Prots[protindex] = ID
                            continue
                        else:
                            break
                    ###ElseIf the accession is found in the Secondary_AC_dict
                    elif Secondary_AC_dict.get(ID, "-1") != "-1":
                        if Prots[protindex] == '':
                            ###Get the corresponding Primary Uniprot Accession ID
                            Prots[protindex] = Secondary_AC_dict[ID]
                            continue
                        else:
                            break
                elif (re_uniprot_missed.match(altID)):
                    sys.exit("ID is a uniprot Accession but failed to grab it for the line:\n" + line)
                elif (re_GeneID.match(altID)):
                    ID = re_GeneID.match(altID).group(1)
                    ##Check if it exists in the dictionary
                    if GeneID_dict.get(ID, False):
                        if Prots[protindex] == '':
                            Prots[protindex] = GeneID_dict[ID]
                            continue
                        else:
                            break
                elif (re_GeneID_missed.match(altID)):
                    sys.exit("ID is a GeneID but failed to grab it for the line:\n" + line)
            if (Prots[0] == '') or (Prots[1] == ''):
                #break_flag = True
                break
        #if break_flag:
            #break
            elif re_IntDetectMethod.match(line_fields[6]):
                IntDetectMethod = re_IntDetectMethod.match(line_fields[6]).group(1)
            elif re_psimi_missed.match(line_fields[6]):
                print("Failed to grab the Interaction Detection Method psi-mi ID for the line:", line)
                break
            elif (re_PMID.match(line_fields[8])):
                PMID = re_PMID.match(line_fields[8]).group(1)
            elif (re_PMID_missed.match(line_fields[8])):
                print("Failed to grab the Pubmed Id for the line", line)
                break
            elif (re_IntType.match(line_fields[11])):
                Interaction_type = re_IntType.match(line_fields[11]).group(1)
            elif (re_psimi_missed.match(line_fields[11])):
                print("Failed to grab the Interaction_type psi-mi Id for the line:", line)
                break
            #else: grabbed all the necessary data, wrtie to output file and move to next line




####Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
    Program: Parses a Protein-protein Interaction (PPI) file, maps to the uniprot file and produces an output file...
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inInteraction', metavar = "Input File", dest = "inInteraction", help = 'Input File Name (Protein-protein Interaction file)', required = True)
    required.add_argument('--inPrimaryAC', metavar = "Input File", dest = "inPrimAC", help = 'Uniprot Primary Accession File generated by the uniprot parser', required = True)
    required.add_argument('--inSecondaryAC', metavar = "Input File", dest = "inSecAC", help = 'Secondary Accession File generated by the uniprot parser', required = True)
    required.add_argument('--inGeneID', metavar = "Input File", dest = "inGeneID", help = 'GeneID File generated by the uniprot parser', required = True)
    #required.add_argument('-o', '--output', metavar = "Output File", dest = "output", help = 'Output File Name', required = True)
    file_parser.set_defaults(func=interaction_parser)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
