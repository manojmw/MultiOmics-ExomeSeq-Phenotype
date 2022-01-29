#!/usr/bin/python

import re
import argparse

###Program to check the number of Human-Human Protein Interactions

def checkHumanPPI(args):

    ##User input -> Protein-protein Interaction file
    interaction_file = open(args.inInteraction)

    ###Skip header
    interaction_file.readline()

    ###Keeping the count of Human-Human Interactions
    HumanInt_Count = 0

    ###Compiling all the regular expressions###

    ###TaxID id for protein
    re_taxID = re.compile('^taxid:(9606)')

    ###Parsing the interaction file
    for line in interaction_file:
        line = line.rstrip('\n')
        line_fields = line.split('\t')

        ###Initializing accumulators
        TaxIDs = ['','']

        for taxindex in [0,1]:
            if (re_taxID.match(line_fields[taxindex+9])):
                ID = re_taxID.match(line_fields[taxindex+9]).group(1)
                TaxIDs[taxindex] = ID
                continue
        if (TaxIDs[0] == '9606') and (TaxIDs[1] == '9606'):
            HumanInt_Count += 1
            continue
    print("\nTotal Number of Human-Human Protein Interactions: ", HumanInt_Count, "\n")        

####Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
-----------------------------------------------------------------------------------------------------------
Program: Parses a MITAB 2.5 or 2.7 file and prints the number of Human-Human Protein Interactions to STDOUT
-----------------------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inInteraction', metavar = "Input File", dest = "inInteraction", help = 'Input File Name (Protein-protein Interaction file)', required = True)
    file_parser.set_defaults(func=checkHumanPPI)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
