## CLI for Chimera_buster
## By: Jessica L Albert
## Last edited : 3/5/24

import argparse
import os
from Chimera_buster.chimera_buster import *

def main():
    # create parser object
    parser = argparse.ArgumentParser(prog = "Chimera_buster",
                                     formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description =('''UMI-based Chimera Buster
Author: Jessica L Albert'''))
 
    # defining arguments for parser object
    parser.add_argument("input_final", type = str, 
                        metavar = "input_final_name", default = None,
                        help = "Designates clusters_concesus.fasta file from the clustering_consensus folder to be filtered. This is required.")

    parser.add_argument("input_prelim", type = str, 
                        metavar = "input_prelim_name", default = None,
                        help = "Designates clusters_concesus.fasta file from the clustering folder to be filtered. This is required.")


    parser.add_argument("output", type = str, 
                        metavar = "output_file_prefix", default = None,
                        help = "Designates output file prefix. This is required.")
    
    parser.add_argument("-m", "--mismatch", type = int,
                        metavar = "int", default = 1,
                        help = "Designates the maximun number of mismatched/indel bases. Default is 1.")

    parser.add_argument("-c", "--check_clusters", type = str,
                        metavar = "True/False", default = "False",
                        help = "When set to True, removed reads whose fwd and rev UMIs both match a kept read are reported as duplicates (the same molecule split by clustering) instead of chimeras. WARNING: this part of the code is slow and is only recommended for low input samples.")
 
    # parse the arguments from standard input
    args = parser.parse_args()

    sample_file = args.input_final

    size_file = args.input_prelim
    
    output_name = args.output
         
         
    # calling functions depending on type of argument
    if args.mismatch !=None:
        mismatch_tolerance = args.mismatch
        
    if args.check_clusters !=None:
        check_clusters_status = "True" if args.check_clusters.lower() == "true" else "False"
             
    chimera_buster(sample_file, size_file, output_name, mismatch_tolerance,check_clusters_status)
    
if __name__ == "__main__":
    # calling the main function
    main()

