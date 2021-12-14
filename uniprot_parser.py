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
    OX = 0
    ENSTs = []
    ENSGs = []
    GeneIDs = []

    ###Compiling all the regular expressions###

    ### Accession Numbers (AC), strip last ';'
    re_AC = re.compile('^AC\s+(\S.*);$')
    ###Organism from the OX line : some lines violate the spec, still grab TAXIDs
    ### even if crap comes after the taxid eg NCBI_TaxID=32201 {ECO:0000312|EMBL:ABW86978.1};
    re_OX = re.compile('^OX\s+NCBI_TaxID=(\d+)[; ]')
    ###DR -> Ensembl transcripts and Genes
    re_ENS = re.compile('^DR\s+Ensembl; (\w+); \w+; (\w+)\.')
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
        elif (re_OX.match(line)):
            if (OX != 0):
                print("Error: Several OX lines for protein\t", ACs)
                break
            OX = re_OX.match(line).group(1)
        elif (re.match(r'^OX\s', line)):
            print("Error: Missed the OX line %s\n", line)
            break
        elif (re_ENS.match(line)):
            ENSTs.append(re_ENS.match(line).group(1))
            ENSGs.append(re_ENS.match(line).group(2))
        elif (re.match(r'^DR\s+Ensembl;', line)): ##If we missed DR->Ensembl lines: die
            print("###Oops....Failed to get all the Ensembl Identifiers\n", ACs, line)
            break
        elif (re_GID.match(line)):
            GeneIDs.append(re_GID.match(line).group(1))
        elif (re.match(r'^DR\s+GeneID.*', line)): ##If any DR line with GeneIDs is missed, we will use this to break the loop
            print("###Oops....Missed the GeneIDs %s\n", GeneID_exists)
            break
        ###Processing the matched records of the protein
        elif (line == '//'):
            # ignore entry if bad species
            if ((OX == 9606) or (OX == 9999)):
                try:
                    ACs_split = ACs.split('; ')
                    primary_AC = ACs_split[0]
                    secondary_ACs = ';'.join(ACs_split[1:])
                    # print("primary:", primary_AC)
                    #  print("secondary:", secondary_ACs)
                    # print to files as needed
                except:
                    print('###Oops...Failed to get Accession IDs')
                    break
                try:
                    # print OX info
                    x=1
                except:
                    print('###Oops...Failed to get the Organism name')
                    break
                try:
                    print(ENSGs)
                except:
                    print('###Oops...Failed to get Ensembl Identifiers')
                    break
                try:
                    GeneIDs_split = [i.split(', |; ') for i in GeneIDs]
                    primary_GeneIDs = str(GeneIDs_split[:1])
                    continue
                except:
                    print('###Oops...Failed to get Ensembl Identifiers')
                    break
            # in any case, reset all accumulators
            ACs = ''
            OX = 0
            ENSTs = []
            ENSGs = []
            GeneIDs = []



except IOError as e:
    print("###Unable to open the file for writing###")
tsv_out.close() ###Closing the written file
