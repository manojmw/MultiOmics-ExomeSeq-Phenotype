#!/usr/bin/python

# manojmw
# 04 Dec, 2021

import re
import argparse
import logging
import sys

###########################################################

# Parses on STDIN a Uniprot file and extracts the following fields from each record:
# - Primary Accession and Secondary Accession(s) from the 'AC' line
# - Gene Name and synonyms from the 'GN' line
# - Taxonomy Identifier from the 'OX' line
# - ENST(s), ENSG(s) & GeneID(s) from the 'DR' line
#
# Prints to STDOUT in .tsv format
# Ouput consists of 7 columns (one record per line):
# - Uniprot Primary Accession
# - Taxonomy Identifier
# - ENST (or a comma seperated list of ENSTs)
# - ENSG (or a comma seperated list of ENSGs)
# - Uniprot Secondary Accession (or a comma seperated list of Uniprot Secondary Accessions)
# - GeneID (or a comma seperated list of GeneIDs)
# - Gene Name (or a comma seperated list of Gene Names)
def uniprot_parser(UniProtinFile):

    logging.info("Starting to run...")

    UniProtinFile = sys.stdin

    # Header line
    UniProt_header = ['Primary_AC', 'TaxID', 'ENSTs', 'ENSGs', 'Secondary_ACs', 'GeneIDs', 'GeneNames']
    print('\t'.join(UniProt_header))

    # Initializing variables/accumulators
    ACs = ''
    TaxID = 0
    ENSTs = []
    ENSGs = []
    GeneIDs = []
    GeneNames = []

    # Compiling all the regular expressions

    # Accession Numbers (AC), strip trailing ';'
    re_AC = re.compile('^AC\s+(\S.*);$')

    # Gene Name/Synonyms from the GN line
    # the below regex pattern ends with either ; or , because 
    # UniProt has a bug (at the time of writing this script) - 
    # the lines for some records are incomplete and does not 
    # end with semi-colon as it is supposed to and even unclosed curly braces
    # Ex: 
    # GN   Name=aadK {ECO:0000303|PubMed:17609790, ECO:0000303|PubMed:8293959,
    # Sometimes, GN lines can also start just with 'Synonyms'
    # Ex:
    # GN   Synonyms=AO {ECO:0000303|PubMed:27255930}
    re_GN = re.compile('^GN\s+((Name=\S.*)|(Synonyms=\S.*))[,;]$')

    # Organism (TaxID) from the OX line; some lines violate the uniprot spec
    # Grab only TaxID even if additional info comes after the TaxID
    # eg NCBI_TaxID=32201 {ECO:0000312|EMBL:ABW86978.1};
    re_TaxID = re.compile('^OX\s+NCBI_TaxID=(\d+)[; ]')

    # Ensembl transcripts and Genes from the DR line
    re_ENS = re.compile('^DR\s+Ensembl; (\S+); \S+; (\S+)\.')

    # GeneIDs from the DR line
    re_GID = re.compile('^DR\s+GeneID;\s+(\d+);')

    # Data lines
    for line in UniProtinFile:
        line = line.rstrip('\r\n') # removing trailing new lines and carriage returns

        # Matching and retrieving the records
        if (re_AC.match(line)):
            if (ACs != ''):
                # Add trailing separator to the previous entries - distingushes each AC
                ACs += '; '
            ACs += re_AC.match(line).group(1)
        elif (re.match(r'^AC\s', line)):
            # If any AC line is missed, Exit the program with an error message
            sys.exit("Error: Missed the AC line %s\n" + line)       
        elif (re_GN.match(line)):
            GNLine = re_GN.match(line).group(1)
            # As per the UniProt documentation, the GN
            # line is supposed to have a defined pattern
            # but that's not the case in reality
            # So we grab everything from the GN line
            # The GN line always ends with a ';' except in some cases where
            # the records are incomplete and has a bug like below (ends with a ',')
            # GN   Name=aadK {ECO:0000303|PubMed:17609790, ECO:0000303|PubMed:8293959,
            try:
                # If the GN line contains any additional info other than
                # Gene Name, it will also be seperated by a '; '
                # Ex: GN   Name=Jon99Cii; Synonyms=SER1, SER5, Ser99Da; ORFNames=CG7877;
                GNList = GNLine.split('; ')
                if GNList:
                    for geneinfo in GNList:
                        # Retrieving only the Gene Name
                        if re.match('^Name=(\S.*)', geneinfo):
                            GeneNameD = re.match('^Name=(\S.*)', geneinfo).group(1)
                            # Sometimes, Gene Name can contain additional info
                            # such as pubmed and other accession IDs
                            # So, we want eliminate these from our
                            # Gene name
                            # Ex: GN   Name=dbaA {ECO:0000303|PubMed:23001671};
                            try:
                                # So we split at the ' {' and keep only the
                                # Gene Name
                                GeneNamewithAcID = GeneNameD.split(' {')
                                GeneName = GeneNamewithAcID[0]
                            except: 
                                GeneName = GeneNameD
                                # There can be multiple 'GN' lines, especially
                                # when there are many synonyms
                                # So we check to avoid adding the same gene name again
                                # Ex: 
                                # GN   Name=Jon99Cii; Synonyms=SER1, SER5, Ser99Da; ORFNames=CG7877;
                                # GN   Name=Jon99Ciii; Synonyms=SER2, SER5, Ser99Db; ORFNames=CG15519;
                            if not GeneName in GeneNames:
                                GeneNames.append(GeneName)   
                        # retreiving synonyms for Gene Name (if it exists)            
                        if re.match('^Synonyms=(\S.*)', geneinfo):
                            GeneSynoinfo = re.match('^Synonyms=(\S.*)', geneinfo).group(1)
                            GeneSynonyms = []
                            # There can be multiple synonyms seperated by a ','
                            # Ex: 
                            # GN   Name=Jon99Cii; Synonyms=SER1, SER5, Ser99Da;
                            try:
                                # Splitting at ',' can sometimes add additional
                                # info of synonym to the list as well
                                # Especially, when more than one additional info
                                # exists for a given Gene synonym and these are also
                                # Seperated by a ','. We eliminate this later while 
                                # adding synonyms to the Gene Names list
                                #
                                # Gene name Synonym with one additional info
                                # Ex: GN   Synonyms=ALI2 {ECO:0000250|UniProtKB:Q80ZD8}
                                #
                                # Gene name Synonym with more than one additional info
                                # seperated by a ','
                                # Ex: GN   Synonyms=OPCL1 {ECO:0000303|PubMed:18267944, ECO:0000303|PubMed:19704801};
                                GeneSynonymList = GeneSynoinfo.split(', ')
                                if GeneSynonymList:
                                    # Like Gene Name, even Gene synonyms can
                                    # can contain additional info such as pubmed and other accession IDs
                                    # So, we want eliminate these and keep only the synonym name
                                    # Ex:
                                    # GN   Name=Sh3bp5; Synonyms=Sab {ECO:0000303|PubMed:10339589};
                                    for synonym in GeneSynonymList:
                                        try:
                                            Gsynowithaddinfo = synonym.split(' {')
                                            Gsynonym = Gsynowithaddinfo[0]
                                            GeneSynonyms.append(Gsynonym)
                                        except:
                                            GeneSynonyms.append(synonym)
                            except:
                                GeneSynonyms.append(GeneSynoinfo)
                            # Avoid adding the same synonym again especially
                            # when they occur on multiple 'GN' lines
                            for synonym in GeneSynonyms:
                                # Eliminating additional info (i.e not a synonym)
                                if not ':' in synonym:
                                    if not synonym in GeneNames:
                                        GeneNames.append(synonym)   
            except:
                # if the GN line contains only the Gene name, we do not split 
                # Ex: GN   Name=APOM; 
                if re.match('^Name=(\S.*)', GNLine):
                    GeneNameD = re.match('^Name=(\S+)', GNLine).group(1)
                    try:
                        # When Gene Name contains additional info
                        GeneNamewithAcID = GeneNameD.split(' {')
                        GeneName = GeneNamewithAcID[0]
                    except:
                        GeneName = GeneNameD 
                    if not GeneName in GeneNames:
                        GeneNames.append(GeneName)     
        elif (re.match(r'^GN\s+Name.*', line)):
            sys.exit("Error: Missed the Gene Name \n" + ACs + line)  
        elif (re.match(r'^GN\s+Synonyms.*', line)):
            sys.exit("Error: Missed the Gene Name Synonym \n" + ACs + line)
        elif (re_TaxID.match(line)):
            if (TaxID != 0):
                sys.exit("Error: Several OX lines for the protein: \t" + ACs)
            TaxID = re_TaxID.match(line).group(1)
        elif (re.match(r'^OX\s',line)):
            sys.exit("Error: Missed the OX line %s\n" + line)
        elif (re_ENS.match(line)):
            ENS_match = re_ENS.match(line)
            ENSTs.append(ENS_match.group(1))
            ENSG = ENS_match.group(2)
            if ENSG not in ENSGs:
                ENSGs.append(ENSG)
        elif (re.match(r'^DR\s+Ensembl;', line)):
            sys.exit("Error: Failed to get all the Ensembl Identifiers\n" + ACs + line)
        elif (re_GID.match(line)):
            GeneIDs.append(re_GID.match(line).group(1))
        elif (re.match(r'^DR\s+GeneID.*', line)):
            sys.exit("Error: Missed the GeneIDs \n" + ACs + line)

        # '//' means End of the record
        # we Process the retreived data
        elif (line == '//'):
            # ignore entry if bad species; TaxID Human = 9606
            if (TaxID == '9606'):
                try:
                    ACs_split = ACs.split('; ')
                    primary_AC = ACs_split[0] # Grab only the first AC
                    secondary_AC_list = ACs_split[1:] # Grab the remaining ACs
                    # storing secondary_ACs as a single comma-seperated string
                    secondary_ACs = ','.join(secondary_AC_list) 
                except:
                    sys.exit('Error: Failed to store Accession IDs for the protein: \t' + ACs)
                try:
                    # storing Gene names as a single comma-seperated string
                    GeneNames = ','.join(GeneNames)
                except:
                    sys.exit('Error: Failed to store Gene Name for the protein: \t' + ACs)    
                try:
                    # storing ENSTs and ENSGs as a single comma-seperated string
                    ENSTs = ','.join(ENSTs)
                    ENSGs = ','.join(ENSGs)
                except:
                    sys.exit('Error: Failed to store Ensembl Identifiers for the protein: \t' + ACs)
                try:
                    # storing GeneIDs as a single comma-seperated string
                    GeneIDs = ','.join(GeneIDs)
                except:
                    sys.exit('Error: Failed to store Gene Identifiers for the protein: \t' + ACs)    

                # Printing to STDOUT
                UniProt_outline = [primary_AC, TaxID, ENSTs, ENSGs, secondary_ACs, GeneIDs, GeneNames]
                print('\t'.join(UniProt_outline))

            # Reset all accumulators and move on to the next record
            ACs = ''
            TaxID = 0
            ENSTs = []
            ENSGs = []
            GeneIDs = []
            GeneNames = []
            continue

    logging.info("All Done, completed successfully!")

    return

###########################################################

# Taking and handling command-line arguments
def main():

    file_parser = argparse.ArgumentParser(description =
    """
--------------------------------------------------------------------------------------------------
Program: Parses on STDIN a UniProt file, processes each record and prints to STDOUT in .tsv format
--------------------------------------------------------------------------------------------------
The output consists of 7 columns:
 -> Uniprot Primary Accession
 -> Taxonomy Identifier
 -> ENST (or a comma seperated list of ENSTs)
 -> ENSG (or a comma seperated list of ENSGs)
 -> Uniprot Secondary Accession (or a comma seperated list of Uniprot Secondary Accessions)
 -> GeneID (or a comma seperated list of GeneIDs)
 -> Gene Name (or a comma seperated list of Gene Names)
--------------------------------------------------------------------------------------------------

Arguments [defaults] -> Can be abbreviated to shortest unambiguous prefixes
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    file_parser.set_defaults(func=uniprot_parser)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    # Logging to Standard Error
    Log_Format = "%(levelname)s - %(asctime)s - %(message)s \n"
    logging.basicConfig(stream = sys.stderr, format  = Log_Format, level = logging.DEBUG)
    main()
