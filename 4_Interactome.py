#!/usr/bin/python

import sys, argparse

##Creating dictionary from the interaction_parser.py output files

def IntPMID(curatedIntfile):

    header = ('Protein_A_UniprotPrimAC', 'Protein_B_UniprotPrimAC', 'Publication_Count', 'Publication_Identifier(s)', 'Experiment_count')
    print('\t'.join(header))

    PPI_PMID_dict = {} ###Dictionary for interactions
    PPI_IntDetMethod_dict = {} ###Dictionary for experiments and filtering based on IntDetMethod

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
            ##Value -> Pubmed Identifier (PMID) - curatedPPI_fields[3]
            (Int_key1, PMID) = (Interactors, curatedPPI_fields[3])

            ##Check if the Key exists in PPI_PMID_dict
            ##If yes, then store the values (PMIDs) as a list
            if PPI_PMID_dict.get(Int_key1, False):
                if PMID not in PPI_PMID_dict[Int_key1]: ###Avoiding duplicate PMIDs
                    PPI_PMID_dict[Int_key1].append(PMID)
            else:
                PPI_PMID_dict[Int_key1] = [PMID]

            ##Key -> UniProt PrimAC of Protein A & B joined together by an '_'
            ##Value -> Interaction Detection Method - curatedPPI_fields[2]
            (Int_key2, IntDetMeth) = (Interactors, curatedPPI_fields[2])

            if PPI_IntDetMethod_dict.get(Int_key2, False):
                PPI_IntDetMethod_dict[Int_key2].append(IntDetMeth)
            else:
                PPI_IntDetMethod_dict[Int_key2] = [IntDetMeth]

    ###Processing the dictionaries and printing to STDOUT
    for (Int_key1,PMID), (Int_key2,IntDetMeth) in zip(PPI_PMID_dict.items(), PPI_IntDetMethod_dict.items()):

        ###Checking if at least one of the experiments for a given interaction has been proved by any binary interaction method
        ###i.e. Other than Affintity chromatography technology (ACT) - 'MI:0004'
        ###Used to eliminate PPIs that have been proved ONLY by using ACT
        if (any(value != 'MI:0004' for value in PPI_IntDetMethod_dict[Int_key2])):
            Proteins = Int_key1.split('_')
            Protein_A = Proteins[0]
            Protein_B = Proteins[1]

            Pubmed_Identifier = ', '.join(PPI_PMID_dict[Int_key1])
            PMID_count = len(PPI_PMID_dict[Int_key1])
            Exp_count = len(PPI_IntDetMethod_dict[Int_key2])

            ###Final Quality Control
            ##At least 2 publications for a given PPI
            ##OR if it is one publication, then at least 3 experiments
            if ((PMID_count >= 2) or (PMID_count == 1 and Exp_count >= 3)):
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

    % cat curatedPPI_file1 curatedPPI_file2 | python 4_Interactome.py

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
