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
        if Secondary_AC_dict.get(key, False):
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
        GeneID_dict[GeneID] = PrimAC
    return GeneID_dict

###Protein-Protein Interaction Parser
def interaction_parser(args):

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

    for line in interaction_file:
        line = line.rstrip('\n')
        line_fields = line.split('\t')

        ###Initializing variables/accumulators
        Prots = ['','']
        Method = ''
        PMID = 0
        Interaction_type = ''

        for protindex in [0,1]:
            if (re_uniprot.match(line_fields[protindex])):
                ID = re_uniprot.match(line_fields[protindex]).group(2)
                ##Check if it exists in the dictionary
                if Primary_AC_dict.get(ID, False):
                    Prots[protindex] = ID
                    continue
            elif (re_uniprot_missed.match(line_fields[protindex])):
                sys.exit("ID is a uniprot Accession but failed to grab it\n" + line)
            elif (re_GeneID.match(line_fields[protindex])):
                ID = re_GeneID.match(line_fields[protindex]).group(1)
                ##Check if it exists in the dictionary
                if GeneID_dict.get(ID, False):
                    Prots[protindex] = GeneID_dict[ID]
                    continue
            elif (re_GeneID_missed.match(line_fields[protindex])):
                sys.exit("ID is a GeneID but failed to grab it\n" + line)

            ###Uniprot AC not found/not primary_AC and GeneID not found, look in Alternate ID columns
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
                            continue next line
                    ###ElseIf the accession is found in the Secondary_AC_dict
                    elif Secondary_AC_dict.get(ID, "-1") != "-1":
                        if Prots[protindex] == '':
                            ###Get the corresponding Primary Uniprot Accession ID
                            Prots[protindex] = Secondary_AC_dict[ID]
                            continue
                        else:
                            continue next line
                elif (re_uniprot_missed.match(altID)):
                    sys.exit("ID is a uniprot Accession but failed to grab it\n" + line)
                elif (re_GeneID.match(altID)):
                    ID = re_GeneID.match(altID).group(1)
                    ##Check if it exists in the dictionary
                    if GeneID_dict.get(ID, False):
                        if Prots[protindex] == '':
                            Prots[protindex] = GeneID_dict[ID]
                            continue
                        else:
                            continue next line
                elif (re_GeneID_missed.match(altID)):
                    sys.exit("ID is a GeneID but failed to grab it\n" + line)






####Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(formatter_class=formatter, description =
    """
    Program: Parses a Protein-protein Interaction (PPI) file, maps to the uniprot file and produces an output file...
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inInteraction',  dest = "inInteraction", help = 'Input File Name (Protein-protein Interaction file)', required = True)
    required.add_argument('--inPrimaryAC',  dest = "inPrimAC", help = 'Uniprot Primary Accession File generated by the uniprot parser', required = True)
    required.add_argument('--inSecondaryAC',  dest = "inSecAC", help = 'Secondary Accession File generated by the uniprot parser', required = True)
    required.add_argument('--inGeneID', dest = "inGeneID", help = 'GeneID File generated by the uniprot parser', required = True)
    #required.add_argument('-o', '--output',  dest = "output", help = 'Output File Name', required = True)
    file_parser.set_defaults(func=interaction_parser)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
