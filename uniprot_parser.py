#!/usr/bin/python

import re
import csv

###Creating output file and opening it for writing ----- ****To be completed later**** --------
outputFile1 = open("uniprot_results_main.tsv", "w")
outputFile2 = open("uniprot_results_secondary.tsv", "w")
try:
    with open("uniprot_results.tsv", 'w', newline = '') as tsv_out:
        csv_writer = csv.writer(tsv_out, delimiter = '\t')
        header = ['AC', 'Organism', 'ENST', 'ENSG', 'GeneID']
        csv_writer.writerow(header)

    ###Opening the file in the read mode
    uniprotfile = open('/path-to-the-file/', 'r')

    ###Initializing variables
    ACs = ''
    Organism = []
    ENSTs = []
    ENSGs = []
    GeneIDs = []

    for line in uniprotfile:
        line = line.rstrip() ##removing any trailing spaces

        ###Matching and retrieving the Accession Number (AC)
        AC_match = re.search(r'^AC\s+(\S.*)', line)
        AC_exists = re.search(r'^AC\s+', line) ##If any AC line is missed, we will use this to break the loop

        ###Matching and retrieving Organism from the OS line
        Organism_match = re.search(r'^OS\s+(.*)\.', line)

        ###Matching and retrieving the Ensemble transcript and Genes
        ENS_match1 = re.search(r'^DR\s+Ensembl;\s+(ENST\d+);\sENSP\d+;\s+(ENSG\d+)\.', line)
        ENS_exists1 = re.search(r'^DR\s+Ensembl;\s+ENST\d+', line) ##If any DR line with ENSTs/ENSGs of the above pattern is missed, we will use this to break the loop

        ###Some records contain only ENSGs with the below regex pattern
        ####We will Match and retrieve lines that contain only Ensembl Genes
        ENS_match2 = re.search(r'^DR\s+Ensembl;\s+(ENSG.*)\.', line)
        ENS_exists2 = re.search(r'^DR\s+Ensembl;\sENSG.*', line) ##If any DR line with ENSGs of the above pattern is missed, we will use this to break the loop

        ###Matching and retrieving the Gene IDs
        GeneID_match = re.search(r'^DR\s+GeneID;\s+(\d+);', line)
        GeneID_exists = re.search(r'^DR\s+GeneID.*', line) ##If any DR line with GeneIDs of the above pattern is missed, we will use this to break the loop

        if AC_match:
            ACs += AC_match.group(1) + ' '
        elif AC_exists:
            print("###Oops....Missed the AC line", AC_exists)
            break
        elif Organism_match:
            Organism.append(Organism_match.group(1))
        elif ENS_match1:
            ENSTs.append(ENS_match1.group(1))
            ENSGs.append(ENS_match1.group(2))
        elif ENS_match2:
            ENSG2 = ENS_match2.group(1)
            if ENSG2 not in ENSGs:
                ENSGs.append(ENSG2)
        elif ENS_exists1 or ENS_exists2:
            print("###Oops....Failed to get all the Ensembl Identifiers", ENS_exists1, ENS_exists2)
            break
        elif GeneID_match:
            GeneIDs.append(GeneID_match.group(1))
        elif GeneID_exists:
            print("###Oops....Missed the GeneIDs", GeneID_exists)
            break

        ###Processing the matched records of the protein
        elif (line == '//'):
            try:
                primary_AC = re.search(r'^(^;\s)+').group(1)
            except:
                print('###Oops...Failed to get Primary Accession')
                break


except IOError as e:
    print("###Unable to open the file for writing###")
tsv_out.close() ###Closing the written file
