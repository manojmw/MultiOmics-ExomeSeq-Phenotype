#!/usr/bin/python

import re
import csv

###Creating output files and opening it for writing
outputFile = open("biogrid.tsv", "w")

try:
    with open("biogrid.tsv", 'w', newline = '') as tsv_out:
        csv_writer = csv.writer(tsv_out, delimiter = '\t')
        header = ['Protein_A', 'Protein_B', 'Method', 'PMID', 'TaxID_A', 'TaxID_B', 'Interaction_type']
        csv_writer.writerow(header)

        ###Opening the file in the read mode
        biogridfile = open(input("Please enter the name of the Biogrid interaction file: "), 'r')












        ###Closing the written files
        tsv1_out.close()
        tsv2_out.close()

except IOError as e:
    print("Error: Unable to open the file for writing")
