#!/usr/bin/python

import re
import argparse, sys

###Function to check the number of Human-Human Protein Interaction experiments###
# Takes miTAB 2.5 or 2.7 file as input
# Keeps the count of experiments where both proteins have TaxID = 9606 i.e Human
# Prints the total no. of Human-Human Protein Interaction experiments to STDOUT
def checkHumanPPI(interaction_file):

    # User input -> Protein-protein Interaction file
    interaction_file = sys.stdin

    # Skip header
    interaction_file.readline()

    # Keeping the count of Human-Human Interaction experiments
    HumanInt_Count = 0

    # Compiling regular expression

    # TaxID id for protein
    re_taxID = re.compile('^taxid:(9606)')

    # Parsing the interaction file
    for line in interaction_file:
        line = line.rstrip('\n')
        line_fields = line.split('\t')

        # Initializing accumulators
        TaxIDs = ['','']

        for taxindex in [0,1]:
            if (re_taxID.match(line_fields[taxindex+9])):
                ID = re_taxID.match(line_fields[taxindex+9]).group(1)
                TaxIDs[taxindex] = ID
                continue
        if (TaxIDs[0] == '9606') and (TaxIDs[1] == '9606'):
            HumanInt_Count += 1
            continue
    print("\nTotal Number of Human-Human Protein Interaction experiments: ", HumanInt_Count, "\n")

    # Closing the file
    interaction_file.close()

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
----------------------------------------------------------------------------------------------------------------------
Program: Parses a miTAB 2.5 or 2.7 file and prints the number of Human-Human Protein Interaction experiments to STDOUT
----------------------------------------------------------------------------------------------------------------------
Usage:

    % python 3_check_HumanPPIExp.py < Input file
                        OR
    % cat Input file | python 3_check_HumanPPIExp.py
----------------------------------------------------------------------------------------------------------------------

    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    file_parser.set_defaults(func=checkHumanPPI)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
