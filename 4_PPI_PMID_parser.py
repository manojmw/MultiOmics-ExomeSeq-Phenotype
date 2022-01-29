#!/usr/bin/python

import argparse

##Creating dictionary from the interaction_parser.py output file

def IntPMID(args):
    header = ('Protein_A_UniprotPrimAC', 'Protein_B_UniprotPrimAC', 'Publication Count', 'Publication_Identifier(s)')
    print('\t'.join(header))

    PPI_PMID_dict = {} ##Initializing an empty dictionary
    curatedPPI_file = open(args.incuratedPPI)

    curatedPPI_file.readline() #Skip header

    ###Parsing the curated interaction file
    for line in curatedPPI_file:
        curatedPPI_fields = line.split('\t')

        #curatedPPI_fields[0] -> Protein_A_UniprotPrimAC
        #curatedPPI_fields[1] -> Protein_B_UniprotPrimAC
        Interactors = curatedPPI_fields[0] + '_' + curatedPPI_fields[1]
        ##Key -> UniProt PrimAC of Protein A & B joined together by an '_'
        ##Value -> Pubmed Identifier (PMID)
        (Int_key, PMID) = (Interactors, curatedPPI_fields[3])
        ##Check if the Key exists in PPI_PMID_dict
        ##If yes, then store the values as a list
        if PPI_PMID_dict.get(Int_key, False):
            PPI_PMID_dict[Int_key].append(PMID)
        else:
            PPI_PMID_dict[Int_key] = [PMID]
    for Int_key in PPI_PMID_dict:

        Proteins = Int_key.split('_')
        Protein_A = Proteins[0]
        Protein_B = Proteins[1]

        Pubmed_Identifier = ', '.join(PPI_PMID_dict[Int_key])
        PMID_count = str(len(PPI_PMID_dict[Int_key]))

        interaction_out_line = (Protein_A, Protein_B, PMID_count, Pubmed_Identifier)
        print('\t'.join(interaction_out_line))
    ###Closing the file
    curatedPPI_file.close()
    return


####Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
--------------------------------------------------------------------------------------------------------
Program: Parses the output file produced by the interaction_parser.py, processes it and prints to STDOUT
--------------------------------------------------------------------------------------------------------
The output consists of four columns in .tsv format:
  -> UniProt Primary Accession of Protein A
  -> UniProt Primary Accession of Protein B
  -> Number of Publications associated with the interaction of the above 2 proteins
  -> PMIDs (or comma seperated list of PMIDs)
--------------------------------------------------------------------------------------------------------
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
