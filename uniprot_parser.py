#!/usr/bin/python

import re
import argparse

###UniProt Parser
def uniprot_parser(args):
    try:
        with open(args.outPrimAC, 'w') as PrimAC_outfile, open(args.outSecAC, 'w') as SecAC_outfile, open(args.outGeneID, 'w') as GeneID_outfile:
            PrimAC_header = ['Primary_AC', 'TaxID', 'ENSTs', 'ENSGs']
            SecAC_header = ['Primary_AC', 'Secondary_ACs']
            GeneID_header = ['Primary_AC', 'Secondary_GeneID']
            print('\t'.join(PrimAC_header), file = PrimAC_outfile)
            print('\t'.join(SecAC_header), file = SecAC_outfile)
            print('\t'.join(GeneID_header), file = GeneID_outfile)

            ###Initializing variables/accumulators
            ACs = ''
            TaxID = 0
            ENSTs = []
            ENSGs = []
            GeneIDs = []

            ###Compiling all the regular expressions###

            ###Accession Numbers (AC), strip trailing ';'
            re_AC = re.compile('^AC\s+(\S.*);$')
            ###Organism (TaxID) from the OX line; some lines violate the uniprot spec
            ###Grab only TaxID even if additional info comes after the TaxID
            ###eg NCBI_TaxID=32201 {ECO:0000312|EMBL:ABW86978.1};
            re_TaxID = re.compile('^OX\s+NCBI_TaxID=(\d+)[; ]')
            ###Ensembl transcripts and Genes from the DR line
            re_ENS = re.compile('^DR\s+Ensembl; (\w+); \w+; (\w+)\.')
            ###GeneIDs from the DR line
            re_GID = re.compile('^DR\s+GeneID;\s+(\d+);')

            ###Open and parse the input file

            for line in open(args.inuniprot):
                line = line.rstrip('\r\n') ##removing trailing new lines and carriage returns

                ###Matching and retrieving the records
                if (re_AC.match(line)):
                    if (ACs != ''):
                        #Add trailing separator to the previous entries - distingushes each AC
                        ACs += '; '
                    ACs += re_AC.match(line).group(1)
                elif (re.match(r'^AC\s', line)): ##If any AC line is missed -> break the loop
                    print("Error: Missed the AC line %s\n", line)
                    break
                elif (re_TaxID.match(line)):
                    if (TaxID != 0):
                        print("Error: Several OX lines for the protein: \t", ACs)
                        break
                    TaxID = re_TaxID.match(line).group(1)
                elif (re.match(r'^OX\s',line)):
                    print("Error: Missed the OX line %s\n", line) ##If any OX line is missed -> break the loop
                    break
                elif (re_ENS.match(line)):
                    ENS_match = re_ENS.match(line)
                    ENSTs.append(ENS_match.group(1))
                    ENSG = ENS_match.group(2)
                    if ENSG not in ENSGs:
                        ENSGs.append(ENSG)
                elif (re.match(r'^DR\s+Ensembl;', line)): ##If any DR line wtih Ensembl IDs is missed -> break the loop
                    print("Error: Failed to get all the Ensembl Identifiers\n", ACs, line)
                    break
                elif (re_GID.match(line)):
                    GeneIDs.append(re_GID.match(line).group(1))
                elif (re.match(r'^DR\s+GeneID.*', line)): ##If any DR line wtih GeneIDs is missed -> break the loop
                    print("Error: Missed the GeneIDs \n", ACs, line)
                    break
                ###Processing the matched records of the protein
                elif (line == '//'):
                    # ignore entry if bad species; Human = 9606, Mouse = 10090
                    if ((TaxID == '9606') or (TaxID == '10090')):
                        try:
                            ACs_split = ACs.split('; ')
                            primary_AC = ACs_split[0] ##Grab only the first AC
                            secondary_ACs = ACs_split[1:] ##Grab the remaining ACs
                        except:
                            print('Error: Failed to store Accession IDs for the protein: \t', ACs)
                            break
                        try:
                            ##Processing ENSTs and ENSGs
                            ENSTs = ','.join(ENSTs)
                            ENSGs = ','.join(ENSGs)
                        except:
                            print('Error: Failed to store Ensembl Identifiers for the protein: \t', ACs)
                            break


                        ###Writing to output files
                        primaryAC_line = [primary_AC, TaxID, ENSTs, ENSGs]
                        print('\t'.join(primaryAC_line), file = PrimAC_outfile)

                        for secondary_AC in secondary_ACs:
                            secondaryAC_line = [primary_AC, secondary_AC]
                            print('\t'.join(secondaryAC_line), file = SecAC_outfile)

                        for GeneID in GeneIDs:
                            GeneID_line = [primary_AC, GeneID]
                            print('\t'.join(GeneID_line), file = GeneID_outfile)

                    #Reset all accumulators and move on to the next record
                    ACs = ''
                    TaxID = 0
                    ENSTs = []
                    ENSGs = []
                    GeneIDs = []
                    continue

        ###Closing the files
        PrimAC_outfile.close()
        SecAC_outfile.close()
        GeneID_outfile.close()

    except IOError as e:
        print("Error: Unable to open the files for writing")

####Taking and handling command-line arguments
def main():

    file_parser = argparse.ArgumentParser(description =
    """
-------------------------------------------------------------------------------------
Program: Parses a uniprot file, processes it and produces the following output files:
-------------------------------------------------------------------------------------
Output File 1:  A tab-seperated file (.tsv) with four columns
                 -> UniProt Primary Accession
                 -> Taxonomy Identifier
                 -> ENST (or Comma seperated list of ENSTs)
                 -> ENSG (or Comma seperated list of ENSGs)

Output File 2:  A tab-seperated file (.tsv) with two columns
                 -> UniProt Secondary Accession
                 -> Corresponding UniProt Primary Accession

Output File 3:  A tab-seperated file (.tsv) with two columns
                 -> GeneID
                 -> Corresponding UniProt Primary Accession
-------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    required = file_parser.add_argument_group('Required arguments')
    optional = file_parser.add_argument_group('Optional arguments')

    required.add_argument('--inuniprot',  metavar = "Input File", dest = "inuniprot", help = 'Input File Name (Uniprot File)', required = True)
    required.add_argument('--outPrimaryAC',  metavar = "Output File", dest = "outPrimAC", help = 'Primary Accession File with ENSTs, ENSGs & TaxID', required = True)
    required.add_argument('--outSecondaryAC', metavar = "Output File", dest = "outSecAC", help = 'Secondary Accession File', required = True)
    required.add_argument('--outGeneID', metavar = "Output File", dest = "outGeneID", help = 'GeneID File', required = True)
    file_parser.set_defaults(func=uniprot_parser)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
