#!/usr/bin/python

import sys, argparse

###Function for building the high-quality INTERACTOME###
# Takes the output files produced by the interaction_parser.py as INPUT
# This file consists of five columns i.e.
# Protein A UniProt PrimAC, Protein B UniProt PrimAC, Interaction Detection Method, Pubmed Identifier and Interaction type
# Processes it by filtering based on Interaction Detection Method and Interaction type
# Final filtering - Each interaction has at least 2 experiments
# Prints to STDOUT
# OUTPUT is a High-quality Human Interactome that consists of five columns in .tsv format
# i.e., UniProt PrimAC, Protein B UniProt PrimAC, Publication_Count, PMID(s) and Experiment_count
def IntPMID(curatedIntfile):

    header = ('Protein_A_UniprotPrimAC', 'Protein_B_UniprotPrimAC', 'Publication_Count', 'Publication_Identifier(s)', 'Experiment_count')
    print('\t'.join(header))

    PPI_PMID_dict = {} # Dictionary for interactions
    PPI_IntDetMethod_dict = {} # Dictionary for experiments and filtering based on IntDetMethod

    curatedIntfile = sys.stdin

    # Parsing the curated interaction file
    for line in curatedIntfile:

        curatedPPI_fields = line.split('\t')

        # Filtering out Interactions based on Interaction Detection Method
        IntDetMethod = curatedPPI_fields[2]
        # MI:0096 -> pull down
        # MI:0254 -> genetic interference
        # MI:0686 -> unspecified method

        ###Filtering out Interactions based on Interaction Type
        IntType = curatedPPI_fields[4].rstrip('\n')
        # MI:0407 -> direct interaction
        # MI:0915 -> physical association

        if IntDetMethod not in ['MI:0096', 'MI:0254', 'MI:0686'] and IntType in ['MI:0407', 'MI:0915']:
            # curatedPPI_fields[0] -> Protein_A_UniprotPrimAC
            # curatedPPI_fields[1] -> Protein_B_UniprotPrimAC
            Interactors = curatedPPI_fields[0] + '_' + curatedPPI_fields[1]

            # Key -> UniProt PrimAC of Protein A & B joined together by an '_'
            # Value -> Pubmed Identifier (PMID) - curatedPPI_fields[3]
            (Int_key, PMIDs) = (Interactors, curatedPPI_fields[3])

            # Check if the Key exists in PPI_PMID_dict
            # If yes, then store the values (PMIDs) as a list
            if PPI_PMID_dict.get(Int_key, False):
                if PMIDs not in PPI_PMID_dict[Int_key]: # Avoiding duplicate PMIDs
                    PPI_PMID_dict[Int_key].append(PMIDs)
            else:
                PPI_PMID_dict[Int_key] = [PMIDs]

            # Key -> UniProt PrimAC of Protein A & B joined together by an '_'
            # Value -> Interaction Detection Method - curatedPPI_fields[2]
            (Int_key, IntDetMeth) = (Interactors, curatedPPI_fields[2])

            if PPI_IntDetMethod_dict.get(Int_key, False):
                PPI_IntDetMethod_dict[Int_key].append(IntDetMeth)
            else:
                PPI_IntDetMethod_dict[Int_key] = [IntDetMeth]

    # Processing the dictionaries and printing to STDOUT
    # Since both the dictionaries have the same keys,
    # We can use the keys from PPI_PMID_dict to iterate through PPI_IntDetMethod_dict
    for Int_key in PPI_PMID_dict:

        # Checking if at least one of the experiments for a given interaction has been proved by any binary interaction method
        # i.e. Other than Affintity Chromatography Technology (ACT) - 'MI:0004'
        # Used to eliminate PPIs that have been proved ONLY by using ACT
        if any(value != 'MI:0004' for value in PPI_IntDetMethod_dict[Int_key]):
            Proteins = Int_key.split('_')
            Protein_A = Proteins[0]
            Protein_B = Proteins[1]

            Pubmed_Identifier = ', '.join(PPI_PMID_dict[Int_key])
            PMID_count = len(PPI_PMID_dict[Int_key])
            Exp_count = len(PPI_IntDetMethod_dict[Int_key])

            # Final Quality Control
            # Each interaction has at least 2 experiments
            if Exp_count >= 2:
                interaction_out_line = (Protein_A, Protein_B, str(PMID_count), Pubmed_Identifier, str(Exp_count))
                print('\t'.join(interaction_out_line))

    return


# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
--------------------------------------------------------------------------------------------------------
Program: Parses the output files produced by the interaction_parser.py, processes it and prints to STDOUT
--------------------------------------------------------------------------------------------------------
Usage:

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

    args = file_parser.parse_args()
    IntPMID(args)

if __name__ == "__main__":
    main()
