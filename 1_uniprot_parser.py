#!/usr/bin/python

# manojmw
# 04 Dec, 2021

import re
import argparse
import logging
import sys

###########################################################

# Parses the Uniprot file and extracts the following fields from each record:
# - Primary Accession and Secondary Accession(s) from the 'AC' line
# - Taxonomy Identifier from the 'OX' line
# - ENST(s), ENSG(s) & GeneID from the 'DR' line
#
# Generates 3 output files:
# Output File 1 -> a tab-seperated file with four columns:
# - UniProt Primary AC
# - TaxID
# - ENST(s)
# - ENSG(s)
#
# Output File 2 -> a tab-seperated file with two columns:
# - UniProt Secondary AC
# - Corresponding UniProt Primary AC
#
# Output File 3 -> a tab-seperated file with two columns:
# - GeneID
# - Corresponding UniProt Primary AC
def uniprot_parser(args):
    try:
        with open(args.outPrimAC, 'w') as PrimAC_outfile, open(args.outSecAC, 'w') as SecAC_outfile, open(args.outGeneID, 'w') as GeneID_outfile:
            PrimAC_header = ['Primary_AC', 'TaxID', 'ENSTs', 'ENSGs']
            SecAC_header = ['Primary_AC', 'Secondary_ACs']
            GeneID_header = ['Primary_AC', 'GeneID']
            print('\t'.join(PrimAC_header), file = PrimAC_outfile)
            print('\t'.join(SecAC_header), file = SecAC_outfile)
            print('\t'.join(GeneID_header), file = GeneID_outfile)

            # Initializing variables/accumulators
            ACs = ''
            TaxID = 0
            ENSTs = []
            ENSGs = []
            GeneIDs = []

            # Compiling all the regular expressions

            # Accession Numbers (AC), strip trailing ';'
            re_AC = re.compile('^AC\s+(\S.*);$')

            # Organism (TaxID) from the OX line; some lines violate the uniprot spec
            # Grab only TaxID even if additional info comes after the TaxID
            # eg NCBI_TaxID=32201 {ECO:0000312|EMBL:ABW86978.1};
            re_TaxID = re.compile('^OX\s+NCBI_TaxID=(\d+)[; ]')

            # Ensembl transcripts and Genes from the DR line
            re_ENS = re.compile('^DR\s+Ensembl; (\w+); \w+; (\w+)\.')

            # GeneIDs from the DR line
            re_GID = re.compile('^DR\s+GeneID;\s+(\d+);')

            logging.info("Processing data from UniProt File")

            logging.info("Writing data to output files...")

            # Data lines
            for line in sys.stdin:
                line = line.rstrip('\r\n') # removing trailing new lines and carriage returns

                # Matching and retrieving the records
                if (re_AC.match(line)):
                    if (ACs != ''):
                        # Add trailing separator to the previous entries - distingushes each AC
                        ACs += '; '
                    ACs += re_AC.match(line).group(1)
                elif (re.match(r'^AC\s', line)):
                    # If any AC line is missed, Exit the program with an error message
                    sys.exit("Error: Missed the AC line %s\n", line)
                elif (re_TaxID.match(line)):
                    if (TaxID != 0):
                        sys.exit("Error: Several OX lines for the protein: \t", ACs)
                    TaxID = re_TaxID.match(line).group(1)
                elif (re.match(r'^OX\s',line)):
                    sys.exit("Error: Missed the OX line %s\n", line)
                elif (re_ENS.match(line)):
                    ENS_match = re_ENS.match(line)
                    ENSTs.append(ENS_match.group(1))
                    ENSG = ENS_match.group(2)
                    if ENSG not in ENSGs:
                        ENSGs.append(ENSG)
                elif (re.match(r'^DR\s+Ensembl;', line)):
                    sys.exit("Error: Failed to get all the Ensembl Identifiers\n", ACs, line)
                elif (re_GID.match(line)):
                    GeneIDs.append(re_GID.match(line).group(1))
                elif (re.match(r'^DR\s+GeneID.*', line)):
                    sys.exit("Error: Missed the GeneIDs \n", ACs, line)

                # '//' means End of the record
                # we Process the retreived data
                elif (line == '//'):
                    # ignore entry if bad species; TaxID Human = 9606, TaxID Mouse = 10090
                    if ((TaxID == '9606') or (TaxID == '10090')):
                        try:
                            ACs_split = ACs.split('; ')
                            primary_AC = ACs_split[0] # Grab only the first AC
                            secondary_ACs = ACs_split[1:] # Grab the remaining ACs
                        except:
                            sys.exit('Error: Failed to store Accession IDs for the protein: \t', ACs)
                        try:
                            # Processing ENSTs and ENSGs
                            ENSTs = ','.join(ENSTs)
                            ENSGs = ','.join(ENSGs)
                        except:
                            sys.exit('Error: Failed to store Ensembl Identifiers for the protein: \t', ACs)

                        # Writing to output files
                        primaryAC_line = [primary_AC, TaxID, ENSTs, ENSGs]
                        print('\t'.join(primaryAC_line), file = PrimAC_outfile)

                        for secondary_AC in secondary_ACs:
                            secondaryAC_line = [primary_AC, secondary_AC]
                            print('\t'.join(secondaryAC_line), file = SecAC_outfile)

                        for GeneID in GeneIDs:
                            GeneID_line = [primary_AC, GeneID]
                            print('\t'.join(GeneID_line), file = GeneID_outfile)

                    # Reset all accumulators and move on to the next record
                    ACs = ''
                    TaxID = 0
                    ENSTs = []
                    ENSGs = []
                    GeneIDs = []
                    continue

        # Closing the files
        PrimAC_outfile.close()
        SecAC_outfile.close()
        GeneID_outfile.close()

    except IOError as e:
        print("Error: Unable to open the files for writing")

    logging.info("Done 🎉")

    return

###########################################################

# Taking and handling command-line arguments
def main():

    file_parser = argparse.ArgumentParser(description =
    """
---------------------------------------------------------------------------------------------
Program: Parses a uniprot file (STDIN), processes it and produces the following output files:
---------------------------------------------------------------------------------------------
Output File 1 (--outPrimaryAC):   A tab-seperated file (.tsv) with four columns
                                   -> UniProt Primary Accession
                                   -> Taxonomy Identifier
                                   -> ENST (or Comma seperated list of ENSTs)
                                   -> ENSG (or Comma seperated list of ENSGs)

Output File 2 (--outSecondaryAC): A tab-seperated file (.tsv) with two columns
                                   -> UniProt Secondary Accession
                                   -> Corresponding UniProt Primary Accession

Output File 3 (--outGeneID):      A tab-seperated file (.tsv) with two columns
                                   -> GeneID
                                   -> Corresponding UniProt Primary Accession
---------------------------------------------------------------------------------------------

Arguments [defaults] -> Can be abbreviated to shortest unambiguous prefixes
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--outPrimAC',  metavar = "Output File", dest = "outPrimAC", help = 'Primary Accession File with ENSTs, ENSGs & TaxID', required = True)
    required.add_argument('--outSecAC', metavar = "Output File", dest = "outSecAC", help = 'Secondary Accession File', required = True)
    required.add_argument('--outGeneID', metavar = "Output File", dest = "outGeneID", help = 'GeneID File', required = True)

    args = file_parser.parse_args()
    uniprot_parser(args)

if __name__ == "__main__":
    # Logging to Standard Error
    Log_Format = "%(levelname)s - %(asctime)s - %(message)s \n"
    logging.basicConfig(stream = sys.stderr, format  = Log_Format, level = logging.DEBUG)
    main()
