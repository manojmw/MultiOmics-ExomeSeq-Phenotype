#!/usr/bin/python

import sys, argparse

###########################################################

# Parses output file produced by the K1 method of MONET tool (DREAM Challenge)
# Extracts clusters with size ≥ 2
#
# Prints to STDOUT
# The output consists of following details for each cluster:
# - Cluster ID and cluster size on the first line
# - Name of the node (one per line)
def ExtrClusters_sizeGT2(K1DREAM_file):

    K1DREAM_file = sys.stdin

    # Cluster counter
    cluster_count = 0

    for line in K1DREAM_file:
        line = line.rstrip('\n') # Remove carriage returns
        line_fields = line.split('\t')

        # Eliminating clusters with size < 2
        if len(line_fields) > 4:
            nodes = line_fields[2:]
            cluster_count += 1
            print("Cluster ID: %d ||" % cluster_count, "Size: %d" % len(nodes) )
            for node in nodes:
                print(node)
            print('\n')
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

The output consists of Cluster ID and size followed by name of the nodes in the cluster
----------------------------------------------------------------------------------------------------------------------------------------------
    """,
    formatter_class = argparse.RawDescriptionHelpFormatter)

    file_parser.set_defaults(func=ExtrClusters_sizeGT2)
    args = file_parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
