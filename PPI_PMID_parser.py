#!/usr/bin/python

import argparse

##Creating dictionary from the interaction_parser output file

def IntPMID(args):
    PPI_PMID_dict = {} ##Initializing an empty dictionary
    curatedPPI_file = open(args.incuratedPPI)

    curatedPPI_file.readline() #Skip header

    for line in curatedPPI_file:
        curatedPPI_fields = line.split('\t')

        #curatedPPI_fields[0] -> Protein_A_UniprotPrimAC
        #curatedPPI_fields[1] -> Protein_B_UniprotPrimAC
        Interactors = curatedPPI_fields[0] + '_' + curatedPPI_fields[1]
        ##Key -> UniProt PrimAC of Protein A&B seperated by an '_'
        ##Value -> Pubmed Identifier (PMID)
        (PPI, PMID) = (Interactors, curatedPPI_fields[3])
        ##Check if the Key exists in PPI_PMID_dict
        ##If yes, then store the values as a list
        if PPI_PMID_dict.get(PPI, False):
            PPI_PMID_dict[PPI].append(PMID)
        else:
            PPI_PMID_dict[PPI] = [PMID]
    for PPI in PPI_PMID_dict:
        #Pubmed_Identifier = ', '.join(PPI_PMID_dict[PPI])
        PMID_count = str(len(PPI_PMID_dict[PPI]))
        interaction_out_line = (PPI, PMID_count)
        print('\t'.join(interaction_out_line))
    return




####Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
--------------------------------------------------------------------------------------------------------
Program: Parses the output file produced by the interaction_parser.py, processes it and prints to STDOUT
--------------------------------------------------------------------------------------------------------
The output consists of two columns in .tsv format:
  -> Protein-Protein Interaction
  -> Number of Pubmed Identifier(s)
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--incuratedPPI', metavar = "Input File", dest = "incuratedPPI", help = 'Input File Name (Curated Protein-protein Interaction file)', required = True)
    file_parser.set_defaults(func=IntPMID)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
