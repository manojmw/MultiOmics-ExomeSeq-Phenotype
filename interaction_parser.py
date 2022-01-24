#!/usr/bin/python

import re
import csv
import argparse
import sys

def interaction_parser(args):
    try:
        with open(args.output, 'w', newline = '') as tsv_out:
            csv_writer = csv.writer(tsv_out, delimiter = '\t')
            header = ['Protein_A', 'Protein_B', 'Method', 'PMID', 'Interaction_type']
            csv_writer.writerow(header)

            ###Creating dictionaries from the uniprot output files###

            ###UniProt Primary Accession dictionary
            Primary_AC_dict = {} ##Initializing an empty dictionary
            Primary_AC_file = open(args.input2)

            Primary_AC_file.readline() ###Skip header

            for line in Primary_AC_file:
                Primary_AC_fields = line.split('\t')
                (key,value) = (Primary_AC_fields[0], Primary_AC_fields[1]) ##Key -> primary accession
                Primary_AC_dict[key] = value

            ###UniProt Secondary Accession dictionary
            Secondary_AC_dict = {} ##Initializing an empty dictionary
            Secondary_AC_file = open(args.input3)

            Secondary_AC_file.readline() ###Skip header

            for line in Secondary_AC_file:
                line = line.rstrip("\n") ##removing carriage returns
                Secondary_AC_fields = line.split('\t')
                (key,value) = (Secondary_AC_fields[1], Secondary_AC_fields[0]) ##Key -> secondary accession
                if Secondary_AC_dict.get(key, False):
                    if Secondary_AC_dict[key] != "-1":
                        Secondary_AC_dict[key]  = "-1"
                    #else: Secondary_AC is already bad => NOOP
                else:
                    Secondary_AC_dict[key] = value

            ###GeneID dictionary
            GeneID_dict = {} ##Initializing an empty dictionary
            GeneID_file = open(args.input4)

            GeneID_file.readline() ###Skip header

            for line in GeneID_file:
                line = line.rstrip("\n") ##removing carriage returns
                GeneID_fields = line.split('\t')
                (key,value) = (GeneID_fields[1], GeneID_fields[0]) ##Key -> GeneID
                GeneID_dict[key] = value

            ###User input -> Protein-protein Interaction file
            interaction_file = open(args.input1)

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
                                    Prots[protindex] = Secondary_AC_dict[ID] ###Get the corresponding Primary Uniprot Accession ID
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












            ###Closing the written files
            tsv_out.close()

    except IOError as e:
            print("Error: Unable to open the file for writing")

####Taking and handling command-line arguments
def main():
    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=70)
    file_parser = argparse.ArgumentParser(formatter_class=formatter, description = "Program: Parses a Protein-protein Interaction (PPI) file, maps to the uniprot file and produces the necessary output files")

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('-i1', '--inInteraction',  dest = "input1", help = 'Input File Name (Protein-protein Interaction file)', required = True)
    required.add_argument('-i2','--inPrimaryAC',  dest = "input2", help = 'Uniprot Primary Accession File generated by the uniprot parser', required = True)
    required.add_argument('-i3', '--inSecondaryAC',  dest = "input3", help = 'Secondary Accession File generated by the uniprot parser', required = True)
    required.add_argument('-i4', '--inGeneID', dest = "input4", help = 'GeneID File generated by the uniprot parser', required = True)
    required.add_argument('-o', '--output',  dest = "output", help = 'Output File Name', required = True)
    file_parser.set_defaults(func=interaction_parser)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
