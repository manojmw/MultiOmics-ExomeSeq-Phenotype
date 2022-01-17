#!/usr/bin/python

import re
import csv
import sys

try:
    with open("biogrid.tsv", 'w', newline = '') as tsv_out:
        csv_writer = csv.writer(tsv_out, delimiter = '\t')
        header = ['Protein_A', 'Protein_B', 'Method', 'PMID', 'TaxID_A', 'TaxID_B', 'Interaction_type']
        csv_writer.writerow(header)

        ###Creating dictionaries from the uniprot output files###

        ###UniProt Primary Accession dictionary
        Primary_AC_dict = {} ##Initializing an empty dictionary
        Primary_AC_file = open(input("Please enter the name of the Uniprot Primary Accession file: "))

        Primary_AC_file.readline() ###Skip header

        for line in Primary_AC_file:
            Primary_AC_fields = line.split('\t')
            (key,value) = (Primary_AC_fields[0], 1) ##Key -> primary accession
            Primary_AC_dict[key] = value

        ###UniProt Secondary Accession dictionary
        Secondary_AC_dict = {} ##Initializing an empty dictionary
        Secondary_AC_file = open(input("Please enter the name of the Uniprot Secondary Accession file: "))

        Secondary_AC_file.readline() ###Skip header

        for line in Secondary_AC_file:
            line = line.rstrip("\n") ##removing carriage returns
            Secondary_AC_fields = line.split('\t')
            (key,value) = (Secondary_AC_fields[1], Secondary_AC_fields[0]) ##Key -> secondary accession
            Secondary_AC_dict[key] = value

        ###UniProt GeneID dictionary
        GeneID_dict = {} ##Initializing an empty dictionary
        GeneID_file = open(input("Please enter the name of the GeneID file: "))

        GeneID_file.readline() ###Skip header

        for line in GeneID_file:
            line = line.rstrip("\n") ##removing carriage returns
            GeneID_fields = line.split('\t')
            (key,value) = (GeneID_fields[1], GeneID_fields[0]) ##Key -> GeneID
            GeneID_dict[key] = value

        ###User input -> Protein-protein Interaction file
        interaction_file = open(input("Please enter the name of the interaction file: "))

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
                elif (Prots[0] == '' or Prots[1] == ''): ###If the protein is not found in the Primary_AC_dict
                    for protindex in [2,3,4,5]:
                        if (re_uniprot.match(line_fields[protindex])):
                            ID = re_uniprot.match(line_fields[protindex]).group(2)
                            ##Check if it exists in the dictionary
                            if ID in Primary_AC_dict or Secondary_AC_dict:
                                Prots[protindex] = ID
                elif (re_uniprot_missed.match(line_fields[protindex])):
                    print("ID is a uniprot Accession but failed to grab it", line)
                    break                


            #try:
            #    Protein_A = re.match('uniprot/swissprot:([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})', line[2]).group(1)

            #except:
            #    print("Failed to retrieve Accession ID for Protein_A")
            #    break













    ###Closing the written files
    tsv_out.close()

except IOError as e:
    print("Error: Unable to open the file for writing")
