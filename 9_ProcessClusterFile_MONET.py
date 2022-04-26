#!/usr/bin/python

# manojmw
# 16 Apr, 2022

import sys, argparse

###########################################################

# Parses output file produced by the MONET tool (DREAM Challenge)
# [Choobdar, Sarvenaz et al. “Assessment of network module identification 
# across complex diseases.” Nature methods vol. 16,9 (2019): 843-852. 
# doi:10.1038/s41592-019-0509-5]
#
# Extracts clusters with size ≥ 2
#
# Prints to STDOUT
# The output consists of following details for each cluster:
# - Cluster ID followed by '||'
# - ENSG of the node (one per line)
#
# End of a given cluster is indicated by an empty line
def ExtrClusters_sizeGT2(DREAM_clusterFile):

    DREAM_clusterFile = sys.stdin

    print("#ClustnSee analysis export")

    # Cluster counter
    cluster_count = 0

    for line in DREAM_clusterFile:
        line = line.rstrip('\n') # Remove carriage returns
        line_fields = line.split('\t')

        # Eliminating clusters with size < 2
        if len(line_fields) >= 4:
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
Program: Parses the output file produced by the MONET tool (DREAM Challenge), extracts clusters with size ≥ 2 and prints to STDOUT
-----------------------------------------------------------------------------------------------------------------------------------------------
Usage:

    % python3 9_ProcessClusterFile_MONET.py < File                         
                        OR    
    % cat File | python3 9_ProcessClusterFile_MONET.py   

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
