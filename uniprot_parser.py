#!/usr/bin/python

import re
import csv

###Creating output file and opening it for reading ----- ****To be completed later**** --------
outputFile1 = open("uniprot_main.tsv", "w")
outputFile2 = open("uniprot_secondary.tsv", "w")

try:
    with open("uniprot_main.tsv", 'w', newline = '') as tsv1_out:
        csv_writer = csv.writer(tsv1_out, delimiter = '\t')
        header = ['Primary_AC', 'TaxID', 'ENST', 'ENSG', 'GeneID']
        csv_writer.writerow(header)

    ###Opening the file in the read mode
    #uniprotfile = open('/Users/macbookpro/Desktop/MBHE/Internship/Data/uniprot_50000.dat', 'r')
    uniprotfile = open(input("Please enter the name of the UniProt file: "), 'r')

    ###Initializing variables/accumulators
    ACs = ''
    TaxID = 0
    ENSTs = ''
    ENSGs = ''
    GeneIDs = []

    ###Compiling all the regular expressions###

    ###Accession Numbers (AC), strip trailing ';'
    re_AC = re.compile('^AC\s+(\S.*);$')
    ###Organism (TaxID) from the OX line; some lines violate the uniprot spec
    ###Grab only TaxID even if additional info comes after the TaxID
    ###eg NCBI_TaxID=32201 {ECO:0000312|EMBL:ABW86978.1};
    re_TaxID = re.compile('^OX\s+NCBI_TaxID=(\d+)[; ]')
    ###Ensembl transcripts and Genes from the DR line
    re_ENS = re.compile('^DR\s+Ensembl; (\w+); \w+; (\w+)\.')
    ###GeneIDs from the DR line
    re_GID = re.compile('^DR\s+GeneID;\s+(\d+);')

    for line in uniprotfile:
        line = line.rstrip('\r\n') ##removing trailing new lines and carriage returns

        ###Matching and retrieving the records
        if (re_AC.match(line)):
            if (ACs != ''):
                #Add trailing separator to the previous entries - distingushes each AC
                ACs += '; '
            ACs += re_AC.match(line).group(1)
        elif (re.match(r'^AC\s', line)): ##If any AC line is missed -> break the loop
            print("Error: Missed the AC line %s\n", line)
            break
        elif (re_TaxID.match(line)):
            if (TaxID != 0):
                print("Error: Several OX lines for the protein: \t", ACs)
                break
            TaxID = re_TaxID.match(line).group(1)
        elif (re.match(r'^OX\s', line)):
            print("Error: Missed the OX line %s\n", line) ##If any OX line is missed -> break the loop
            break
        elif (re_ENS.match(line)):
            ENSTs.append(re_ENS.match(line).group(1))
            ENSGs.append(re_ENS.match(line).group(2))
        elif (re.match(r'^DR\s+Ensembl;', line)): ##If any DR line wtih Ensembl IDs is missed -> break the loop
            print("Error: Failed to get all the Ensembl Identifiers\n", ACs, line)
            break
        elif (re_GID.match(line)):
            GeneIDs.append(re_GID.match(line).group(1))
        elif (re.match(r'^DR\s+GeneID.*', line)): ##If any DR line wtih GeneIDs is missed -> break the loop
            print("Error: Missed the GeneIDs \n", ACs, line)
            break
        ###Processing the matched records of the protein
        elif (line == '//'):
            # ignore entry if bad species
            if ((TaxID == '9606') or (TaxID == '10090')):
                try:
                    ACs_split = ACs.split('; ')
                    primary_AC = ACs_split[0] ##Grab only the first AC
                    secondary_ACs = ';'.join(ACs_split[1:]) ##Grab the remaining ACs
                except:
                    print('Error: Failed to store Accession IDs for the protein: \t', ACs)
                    break
                try:
                    ##Processing ENSTs and ENSGs
                    ENSTs = ','.join(ENSTs)
                    ENSGs = ','.join(ENSGs)
                except:
                    print('Error: Failed to store Ensembl Identifiers for the protein: \t', ACs)
                    break
                try:
                    GeneIDs = ','.join(GeneIDs)
                    print(GeneIDs)
                except:
                    print('Error: Failed to store GeneID for the protein: \t', ACs)
                    break
                ###Writing to the file
                newline1 = [primary_AC, TaxID, ENSTs, ENSGs, GeneIDs]
                newline2 = [primary_AC, secondary_ACs]
                csv_writer.writerow(newline1)
                #csv_writer.writerow(newline2)
            #Reset all accumulators and move on to the next record
            ACs = ''
            TaxID = 0
            ENSTs = []
            ENSGs = []
            GeneIDs = []
            continue

except IOError as e:
    print("Error: Unable to open the file for writing")
tsv1_out.close() ###Closing the written file
