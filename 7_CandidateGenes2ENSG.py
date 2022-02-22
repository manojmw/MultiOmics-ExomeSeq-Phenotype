#!/usr/bin/python

import argparse
import pandas as pd

###Function for extracting data from patient metadata file(s)###
# Takes patient metadata file(s) in .xlsx format as INPUT
# Extracts data from 3 columns - Gene, pathologyID and Confidence score
# Prints to STDOUT in .tsv format
def metaParser(candidate_files):

    # Initializing an empty data frame to store the data from all files
    meta_data = pd.DataFrame()

    # Iterating over the list of files and appending to the DataFrame (meta_data)
    for file in candidate_files:
        data = pd.read_excel(file)
        meta_data = pd.concat([meta_data, data])

    # Extract Gene, pathologyID and Confidence score and drop rows with missing values
    Candidate_data = pd.DataFrame(meta_data, columns=['Gene', 'pathologyID', 'Confidence score']).dropna()
    print(Candidate_data.to_csv(sep = '\t', index = False))

    return


####Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
------------------------------------------------------------------------------
Program: Parses the patient metadata file(s), processes it and prints to STDOUT
------------------------------------------------------------------------------
The output consists of 3 columns in .tsv format:
 -> Gene name
 -> Pathology Identifier
 -> Confidence score
-----------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inMetaFile', metavar = "Input File", dest = "inMetaFile", nargs = '+', help = 'Input File Name (Patient meta data file in xlsx format)', required = True)

    args = file_parser.parse_args()
    metaParser(args.inMetaFile)

if __name__ == "__main__":
    main()
