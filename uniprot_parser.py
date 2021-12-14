#!/usr/bin/python

import re
import csv

###Creating output file and opening it for reading ----- ****To be completed later**** --------
outputFile1 = open("uniprot_results_main.tsv", "w")
outputFile2 = open("uniprot_results_secondary.tsv", "w")

try:
    with open("uniprot_results.tsv", 'w', newline = '') as tsv_out:
        csv_writer = csv.writer(tsv_out, delimiter = '\t')
        header = ['AC', 'Organism', 'ENST', 'ENSG', 'GeneID']
        csv_writer.writerow(header)

    ###Opening the file in the read mode
    #uniprotfile = open('/Users/macbookpro/Desktop/MBHE/Internship/Data/uniprot_50000.dat', 'r')
    uniprotfile = open('sample_uniprotfile.dat', 'r')

    ###Initializing variables
    ACs = ''
    OS = ''
    ENSTs = []
    ENSGs = []
    GeneIDs = []

    ###Compiling all the regular expressions###

    ### Accession Numbers (AC), strip last ';'
    re_AC = re.compile('^AC\s+(\S.*);$')
    ###Compliling Regex for the Organism from the OS line
    re_OS = re.compile('^OS\s+(.*)')
    ###DR -> Ensembl transcripts and Genes
    re_ENS1 = re.compile('^DR\s+Ensembl;\s+(ENST\d+);\sENSP\d+;\s+(ENSG\d+)\.')
    ###Some records contain only ENSGs with the below regex pattern
    ####We will Compile Regex for the lines that contain only Ensembl Genes from the DR line
    re_ENS2 = re.compile('^DR\s+Ensembl;\s+(ENSG.*)\.')
    ###Compliling Regex for the the Gene IDs from the DR line
    re_GID = re.compile('^DR\s+GeneID;\s+(\d+);')

    for line in uniprotfile:
        line = line.rstrip('\r\n') ##removing new lines and carriage returns

        ###Matching and retrieving the records
        if (re_AC.match(line)):
            if (ACs != ''):
                # add trailing separator to previous entries
                ACs += '; '
            ACs += re_AC.match(line).group(1)
        elif (re.match(r'^AC\s', line)): ##If any AC line is missed, we will use this to break the loop
            print("###Oops....Missed the AC line %s\n", line)
            break
        elif (re_OS.match(line)):
            OS += re_OS.match(line).group(1) + '; '
        elif (re.match(r'^OS\s+', line)):
            print("###Oops....Missed the OS line %s\n", line)
            break
        elif (re_ENS1.match(line)):
            ENSTs.append(re_ENS1.match(line).group(1))
            ENSGs.append(re_ENS1.match(line).group(2))
        elif (re_ENS2.match(line)):
            ENSG2 = re_ENS2.match(line).group(1)
            if ENSG2 not in ENSGs:
                ENSGs.append(ENSG2)
        elif (re.match(r'^DR\s+Ensembl;\s+(ENSG.*)\.', line) or re.match(r'^DR\s+Ensembl;\sENSG.*', line)): ##If any DR line with ENSTs/ENSGs is missed, we will use this to break the loop
            print("###Oops....Failed to get all the Ensembl Identifiers %s\n", ENS_exists1, ENS_exists2)
            break
        elif (re_GID.match(line)):
            GeneIDs.append(re_GID.match(line).group(1))
        elif (re.match(r'^DR\s+GeneID.*', line)): ##If any DR line with GeneIDs is missed, we will use this to break the loop
            print("###Oops....Missed the GeneIDs %s\n", GeneID_exists)
            break
        ###Processing the matched records of the protein
        elif (line == '//'):
            try:
                ACs_split = ACs.split('; ')
                primary_AC = ACs_split[0]
                secondary_ACs = ACs_split[1:]
                # print to files as needed
                ACs = ''
            except:
                print('###Oops...Failed to get Accession IDs')
                break
            try:
                OS_split = OS.split(';')
                Organism = str(OS_split[:1])
                OS = ''
            except:
                print('###Oops...Failed to get the Organism name')
                break
            try:
                ENSTs_split = [i.split(', |; ') for i in ENSTs]
                primary_ENST = str(ENSTs_split[:1])
                alternate_ENST = str(ENSTs_split[1:-1])
                ENSGs_split = [i.split(', |; ') for i in ENSGs]
                primary_ENSG = str(ENSGs_split[:1])
                alternate_ENSG = str(ENSGs_split[1:-1])
                ENSTs = []
                ENSGs = []
            except:
                print('###Oops...Failed to get Ensembl Identifiers')
                break
            try:
                GeneIDs_split = [i.split(', |; ') for i in GeneIDs]
                primary_GeneIDs = str(GeneIDs_split[:1])
                GeneIDs = []
                continue
            except:
                print('###Oops...Failed to get Ensembl Identifiers')
                break




except IOError as e:
    print("###Unable to open the file for writing###")
tsv_out.close() ###Closing the written file
