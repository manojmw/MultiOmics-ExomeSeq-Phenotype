#!/usr/bin/python

# manojmw
# 16 Apr, 2022

import sys, argparse

###########################################################

# Parses output file produced by the K1 method of MONET tool (DREAM Challenge)
# Extracts clusters with size ≥ 2
#
# Prints to STDOUT
# The output consists of following details for each cluster:
# - Cluster ID followed by '||'
# - ENSG of the node (one per line)
#
# End of a given cluster is indicated by an empty line
def ExtrClusters_sizeGT2(K1DREAM_file):

    K1DREAM_file = sys.stdin

    print("#ClustnSee analysis export")

    # Cluster counter
    cluster_count = 0

    for line in K1DREAM_file:
        line = line.rstrip('\n') # Remove carriage returns
        line_fields = line.split('\t')

        # Eliminating clusters with size < 2
        if len(line_fields) >= 3:
            nodes = line_fields[2:]
            cluster_count += 1
            print("ClusterID:%d||" % cluster_count)
            for node in nodes:
                print(node)
            print("") 
    return

###########################################################

# Taking and handling command-line arguments
def main():
    file_parser = argparse.ArgumentParser(description =
    """
-----------------------------------------------------------------------------------------------------------------------------------------------
Program: Parses the output file produced by the K1 method of MONET tool (DREAM Challenge), extracts clusters with size ≥ 2 and prints to STDOUT
-----------------------------------------------------------------------------------------------------------------------------------------------
Usage:

    % python 6_K1DREAM_ClustersSizeGT2.py < Input file
                        OR
    % cat Input file | python 6_K1DREAM_ClustersSizeGT2.py

The output consists of:
 -> Cluster ID followed by '||'
 -> ENSG of the node (one per line)
 -> End of a given cluster is indicated by an empty line
----------------------------------------------------------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    file_parser.set_defaults(func=ExtrClusters_sizeGT2)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
