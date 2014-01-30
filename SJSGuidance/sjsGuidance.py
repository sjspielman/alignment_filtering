import argparse
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser(prefix_chars='+-', usage='--protein_file <User File>')
    parser.add_argument("-infile", help="A file containing unaligned AA sequences in FASTA format", required=False, dest="infile", type=str)
    ##parser.add_argument("-format", help="Specifies the format of the input file (fasta or nex", dest="form", type=str, default="fasta")
    parser.add_argument("-n", type=int, help="Number of processes to use", default=1)
    parser.add_argument("-bootstraps", help="The number of bootstraps to perform", required=False,
            dest="bootstraps", default=10)
    parser.add_argument("-alphabet", help="Whether AAs or NTs are used", type=str,
            default="AA", required=False) ##AA or NT, default is AA
    ## Gap penalization is now hard-coded by default in accordance with the original runs
    ## parser.add_argument("-gaps", help="Type of gap penalization", default=0,
    ##        type=int, dest="pflag")

    return parser.parse_args()

def main():
    args = parse_args()
    user_file = args.infile
    while args.infile is None:
        args.infile = raw_input("Please provide a protein file in FASTA or PHYLIP format: ")
    while args.form is None:
        args.form = 
    

main()
