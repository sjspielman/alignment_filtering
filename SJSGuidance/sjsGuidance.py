import subprocess
import argparse
import sys
import os
try:
    from dendropy import *
except:
    print "Please install Dendropy. See the README for details."

def parse_args():
    parser = argparse.ArgumentParser(prefix_chars='+-', usage='--protein_file <User File>')
    parser.add_argument("-infile",help="A file containing unaligned AA sequences in FASTA format", dest="infile", type=str, required=False)
    parser.add_argument("-p",     help="Number of processes to use",                               dest="threads", type=int, required=False, default=1)
    parser.add_argument("-form",  help="The file format (usually FASTA)",                          dest="form", type=str, default="FASTA")
    parser.add_argument("-n",     help="The number of bootstraps to perform",                      dest="bootstraps", type=int, required=False, default=100)
    parser.add_argument("-d",     help="Whether protein or nucleotides are used (prot or dna)",    dest="alphabet", type=str, required=False, default="prot")

    return parser.parse_args()

def getMafft():
    found = subprocess.call(["which", "mafft"])
    if not found == 0:
        print "MAFFT needs to be installed and on the system path"
        print "See the README on how to do this."

def getFastTree():
    found = subprocess.call(["which", "fasttree"])
    if not found == 0:
        print "FastTree needs to be installed and on the system path"
        print "See the README on how to do this"

def getRAxML():
    found = subprocess.call(["which", "raxmlHPC"])
    if not found:
        print "RAxML needs to be installed and on the system path"
        print "See the README on how to do this."
        

def main():
    getMafft()
    getFastTree()
    getRAxML()
    args = parse_args()
    user_file = args.infile
    while args.infile is None:
        args.infile = raw_input("Please provide an input sequence file in FASTA format: ")
    while args.form is None:
        print ""
        args.form = raw_input("Please tell me what format the infile is in.\nIt should be a FASTA file: ")
    if args.threads is None:
        print ""
        print "One thread will be used. To change the number of threads, use the -p flag"
        print "More threads will run faster, but you shouldn't use more than the number of cores in your machine\n"
    if args.bootstraps is None:
        print ""
        print "100 bootstraps will be performed. To change the number of bootstraps, use the -n flag\n"
    if args.alphabet is None:
    	print ""
        print "No alphabet was selected, so amino acids will be used by default. Use the -alphabet flag to specify an alphabet.\n"
    command = "python main.py " + str(args.infile) + " " + str(args.alphabet) + " " + str(args.bootstraps) + " " + str(args.threads)
    subprocess.call(command, shell=True)
    
main()
