#!/usr/bin/python

import sys, argparse

##Creating dictionary from the interaction_parser.py output file

def IntPMID(curatedIntfile):

    header = ('Protein_A_UniprotPrimAC', 'Protein_B_UniprotPrimAC', 'Publication_Count', 'Publication_Identifier(s)', 'Experiment_count')
    print('\t'.join(header))

    PPI_PMID_dict = {} ###Dcitionary for interactions
    PPI_Exp_dict = {} ###Dcitionary for experiments

    curatedIntfile = sys.stdin

    ###Parsing the curated interaction file
    for line in curatedIntfile:

        curatedPPI_fields = line.split('\t')

        ###Filtering out Interactions based on Interaction Detection Method
        IntDetMethod = curatedPPI_fields[2]
        ##MI:0096 -> pull down
        ##MI:0254 -> genetic interference
        ##MI:0686 -> unspecified method

        ###Filtering out Interactions based on Interaction Type
        IntType = curatedPPI_fields[4].rstrip('\n')
        ##MI:0407 -> direct interaction
        ##MI:0915 -> physical association

        if IntDetMethod not in ['MI:0096', 'MI:0254', 'MI:0686'] and IntType in ['MI:0407', 'MI:0915']:
            #curatedPPI_fields[0] -> Protein_A_UniprotPrimAC
            #curatedPPI_fields[1] -> Protein_B_UniprotPrimAC
            Interactors = curatedPPI_fields[0] + '_' + curatedPPI_fields[1]

            ##Key -> UniProt PrimAC of Protein A & B joined together by an '_'
            ##Value -> Pubmed Identifier (PMID)
            (Int_key, PMID) = (Interactors, curatedPPI_fields[3])

            ##Check if the Key exists in PPI_PMID_dict
            ##If yes, then store the values (PMIDs) as a list
            if PPI_PMID_dict.get(Int_key, False):
                if PMID not in PPI_PMID_dict[Int_key]: ###Avoiding duplicate PMIDs
                    PPI_PMID_dict[Int_key].append(PMID)
            else:
                PPI_PMID_dict[Int_key] = [PMID]

            if PPI_Exp_dict.get(Int_key, False):
                PPI_Exp_dict[Int_key].append(PMID) ##Not Avoiding duplicate PMIDs to keep the count of experiments
            else:
                PPI_Exp_dict[Int_key] = [PMID]

    ###Processing the dictionary and printing to STDOUT
    for Int_key, Int_key in zip(PPI_PMID_dict, PPI_Exp_dict):

        Proteins = Int_key.split('_')
        Protein_A = Proteins[0]
        Protein_B = Proteins[1]

        Pubmed_Identifier = ', '.join(PPI_PMID_dict[Int_key])
        PMID_count = len(PPI_PMID_dict[Int_key])
        Exp_count = len(PPI_Exp_dict[Int_key])

        if (PMID_count >= 2) or (PMID_count == 1 and Exp_count >= 3):
            interaction_out_line = (Protein_A, Protein_B, str(PMID_count), Pubmed_Identifier, str(Exp_count))
            print('\t'.join(interaction_out_line))

    return


####Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
--------------------------------------------------------------------------------------------------------
Program: Parses the output files produced by the interaction_parser.py, processes it and prints to STDOUT
--------------------------------------------------------------------------------------------------------
The script should be run using the below command:

    % cat curatedPPI_file1 curatedPPI_file2 | python 4_PPI_PMID_parser.py

-> curatedPPI_files: output files produced by the interaction_parser.py

The output (High-quality Human Interactome) consists of five columns in .tsv format:
  -> UniProt Primary Accession of Protein A
  -> UniProt Primary Accession of Protein B
  -> Number of Publications associated with the interaction of the above 2 proteins
  -> PMID (or a comma seperated list of PMIDs)
  -> Count of Experiments for each interaction
--------------------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    file_parser.set_defaults(func=IntPMID)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
