#!/usr/bin/python

import re
import csv
import argparse

def interaction_parser(args):
    try:
        with open("biogrid.tsv", 'w', newline = '') as tsv_out:
            csv_writer = csv.writer(tsv_out, delimiter = '\t')
            header = ['Protein_A', 'Protein_B', 'Method', 'PMID', 'TaxID_A', 'TaxID_B', 'Interaction_type']
            csv_writer.writerow(header)

            ###Creating dictionaries from the uniprot output files###

            ###UniProt Primary Accession dictionary
            Primary_AC_dict = {} ##Initializing an empty dictionary
            Primary_AC_file = open(args.input2)

            Primary_AC_file.readline() ###Skip header

            for line in Primary_AC_file:
                Primary_AC_fields = line.split('\t')
                (key,value) = (Primary_AC_fields[0], 1) ##Key -> primary accession
                Primary_AC_dict[key] = value

            ###UniProt Secondary Accession dictionary
            Secondary_AC_dict = {} ##Initializing an empty dictionary
            Secondary_AC_file = open(args.input3)

            Secondary_AC_file.readline() ###Skip header

            for line in Secondary_AC_file:
                line = line.rstrip("\n") ##removing carriage returns
                Secondary_AC_fields = line.split('\t')
                (key,value) = (Secondary_AC_fields[1], Secondary_AC_fields[0]) ##Key -> secondary accession
                Secondary_AC_dict[key] = value

            ###UniProt GeneID dictionary
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

            ###Initializing variables/accumulators
            Prots = ['','']
            Method = ''
            PMID = 0
            TaxID_A = 0
            TaxID_B = 0
            Interaction_type = ''

            ###Compiling all the regular expressions###

            ###uniprot ids for protein
            re_uniprot = re.compile('^uniprot(kb|/swiss-prot):([A-Z0-9]+)$')
            re_uniprot_missed = re.compile('^uniprot')
            re_GeneID = re.compile('^entrez gene/locuslink:(\d+)$')
            re_GeneID_missed = re.compile('^entrez gene/locuslink')

            for line in interaction_file:
                line = line.rstrip('\n')
                line_fields = line.split('\t')

                for protindex in [0,1]:
                    if (re_uniprot.match(line_fields[protindex])):
                        ID = re_uniprot.match(line_fields[protindex]).group(2)
                        ##Check if it exists in the dictionary
                        if ID in Primary_AC_dict:
                            Prots[protindex] = ID
                    elif (re_uniprot_missed.match(line_fields[protindex])):
                        print("ID is a uniprot Accession but failed to grab it", line)
                        break
                if (Prots[0] == ''): ###If the Protein_A Accession ID is not found, search in the 3rd column
                    line_fields_split = line_fields[2].split('|')
                    if (re_uniprot.match(line_fields_split)):
                        ID = re_uniprot.match(line_fields_split).group(2)
                        ##Check if it exists in the dictionary
                        if ID in Primary_AC_dict:
                            Prots[0] = ID
                        ###ElseIf the accession is found in the Secondary_AC_dict
                        elif ID in Secondary_AC_dict:
                            Prots[0] = Secondary_AC_dict.get(ID) ###Get the corresponding Primary Uniprot Accession ID
                    elif (re_uniprot_missed.match(line_fields_split)):
                        print("ID is a uniprot Accession but failed to grab it", line)
                        break
                elif (Prots[1] == ''): ###If the Protein_B Accession ID is not found, search in the 4th column
                    line_fields_split = line_fields[3].split('|')
                    if (re_uniprot.match(line_fields_split)):
                        ID = re_uniprot.match(line_fields_split).group(2)
                        ##Check if it exists in the dictionary
                        if ID in Primary_AC_dict:
                            Prots[1] = ID
                        ###ElseIf the accession is found in the Secondary_AC_dict
                        elif ID in Secondary_AC_dict:
                            Prots[1] = Secondary_AC_dict.get(ID) ###Get the corresponding Primary Accession
                    elif (re_uniprot_missed.match(line_fields_split)):
                        print("ID is a uniprot Accession but failed to grab it", line)
                        break
                elif (Prots[0] == ''): ###If the Protein_A Accession ID is not found, search in the 5th column
                    line_fields_split = line_fields[4].split('|')
                    if (re_uniprot.match(line_fields_split)):
                        ID = re_uniprot.match(line_fields_split).group(2)
                        ##Check if it exists in the dictionary
                        if ID in Primary_AC_dict:
                            Prots[0] = ID
                        ###ElseIf the accession is found in the Secondary_AC_dict
                        elif ID in Secondary_AC_dict:
                            Prots[0] = Secondary_AC_dict.get(ID) ###Get the corresponding Primary Uniprot Accession ID
                    elif (re_uniprot_missed.match(line_fields_split)):
                        print("ID is a uniprot Accession but failed to grab it", line)
                        break
                elif (Prots[1] == ''): ###If the Protein_B Accession ID is not found, search in the 5th column
                    line_fields_split = line_fields[5].split('|')
                    if (re_uniprot.match(line_fields_split)):
                        ID = re_uniprot.match(line_fields_split).group(2)
                        ##Check if it exists in the dictionary
                        if ID in Primary_AC_dict:
                            Prots[1] = ID
                        ###ElseIf the accession is found in the Secondary_AC_dict
                        elif ID in Secondary_AC_dict:
                            Prots[1] = Secondary_AC_dict.get(ID) ###Get the corresponding Primary Accession
                    elif (re_uniprot_missed.match(line_fields_split)):
                        print("ID is a uniprot Accession but failed to grab it", line)
                        break
                elif (Prots[0] == ''): ###If the Protein_A Uniprot Accession ID is not found, search using GeneID
                    if (re_GeneID.match(line_fields[0])):
                        ID = re_GeneID.match(line_fields[0]).group(1)
                        ##Check if it exists in the GeneID_dict
                        if ID in GeneID_dict:
                            Prots[0] = GeneID_dict.get(ID) ###Get the corresponding Primary Uniprot Accession ID
                    elif (re_GeneID_missed.match(line_fields[0])):
                        print("Failed to grab GeneID for the line: ", line)
                        break
                elif (Prots[1] == ''): ###If the Protein_A Uniprot Accession ID is not found, search using GeneID
                    if (re_GeneID.match(line_fields[1])):
                        ID = re_GeneID.match(line_fields[1]).group(1)
                        ##Check if it exists in the GeneID_dict
                        if ID in GeneID_dict:
                            Prots[1] = GeneID_dict.get(ID) ###Get the corresponding Primary Uniprot Accession ID
                    elif (re_GeneID_missed.match(line_fields[1])):
                        print("Failed to grab GeneID for the line: ", line)
                        break











        ###Closing the written files
        tsv_out.close()

    except IOError as e:
        print("Error: Unable to open the file for writing")

####Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description = "Program: Parses a Protein-protein Interaction (PPI) file, maps to the uniprot file and produces the necessary output files")
    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=70)
    file_parser = argparse.ArgumentParser(formatter_class=formatter)

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
