#!/usr/bin/python

import sys, argparse

def PPINetworkGen(highqualityPPI):

    highqualityPPI = sys.stdin

    highqualityPPI.readline() #Skip header

    ###Parsing the file
    for line in highqualityPPI:
        line = line.rstrip('\n') ##remove carriage returns
        Interactome_fields = line.split('\t')

        ###Defining a formula that can be used to assign a score for each Interaction
        ###This considers the no. of Publications & Experiment count
        ###The score will be used as the edge weight

        ###Interactome_fields[4] -> Experiment count
        ###Interactome_fields[2] -> PMID count
        ###For every 3 experiments, we assign a score = 1
        ###Each PMID gets a score = 1
        Edge_attrib = int(Interactome_fields[2]) + int(Interactome_fields[4])/3

        ###Printing to STDOUT
        ###Interactome_fields[0] -> Protein A Uniprot Primary Accession
        ###Interactome_fields[1] -> Protein B Uniprot Primary Accession

        print(Interactome_fields[0], "\t", Interactome_fields[1], "\t", round(Edge_attrib, 2))

    return


####Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
-----------------------------------------------------------------------------------------------------------
Program: Parses the output file produced by Interactome.py, processes it and prints to STDOUT
-----------------------------------------------------------------------------------------------------------
The script should be run using one of the below commands:

    % python 5_ModuleInputFile.py < Input file
                        OR
    % cat Input file | python 5_ModuleInputFile.py

You can also direct the output of Interactome.py to this script using the below command:

    % cat curatedPPI_file1 curatedPPI_file2 | python 4_Interactome.py | python 5_ModuleInputFile.py

-> curatedPPI_files: output files produced by the interaction_parser.py

The output consists of three columns in .tsv format:
  -> UniProt Primary Accession of Protein A
  -> UniProt Primary Accession of Protein B
  -> Edge attribute (weight of the edge based on the scores)

The OUTPUT FILE generated by this script can be used as INPUT for most of the MODULE IDENTIFICATION METHODS
-----------------------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    file_parser.set_defaults(func=PPINetworkGen)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
