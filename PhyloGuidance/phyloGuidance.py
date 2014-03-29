import subprocess
import argparse
import sys
import os
try:
    from dendropy import *
except:
    print "Please install Dendropy. See the README for details."

def parse_args():
    parser = argparse.ArgumentParser(prefix_chars='+-', usage='Enter python phyloGuidance.py --help/-h for instructions')
    parser.add_argument("-infile", help="A file containing unaligned sequences in FASTA format", required=False, dest="infile", type=str)
    parser.add_argument("-n",dest="threads", type=int, help="Number of processes to use")
    parser.add_argument("-form",dest="form", type=str, help="The infile format (usually FASTA)", default="FASTA")
    parser.add_argument("-bootstraps", help="The number of bootstraps to perform", required=False,
            dest="bootstraps")
    parser.add_argument("-alphabet", help="Whether AAs or NTs are used (protein or nucleotide)", type=str,
            default="protein", required=False, dest="alphabet") ##AA or NT, default is AA
    ## Gap penalization is now hard-coded by default in accordance with the original runs

    return parser.parse_args()

def getMafft():
    found = subprocess.call(["which", "mafft"])
    assert (found == 0), "MAFFT needs to be installed and on the system path. See the README for instructions."

def getFastTree():
    found = subprocess.call(["which", "fasttree"])
    assert (found == 0), "FastTree needs to be installed and on the system path. See the README for instructions."

def getRAxML():
    found = subprocess.call(["which", "raxmlHPC"])
    assert (found == 0), "RAxML needs to be installed and on the system path. See the README for instructions."

def main():
    args = parse_args()
    getMafft()
    getFastTree()
    getRAxML()
    user_file = args.infile
    while args.infile is None:
        args.infile = raw_input("Please provide a protein file in FASTA format: ")
    while args.form is None:
        print ""
        args.form = raw_input("Please tell me what format the infile is in.\nIt should be a FASTA file: ")
    if args.threads is None:
        print ""
        import multiprocessing 
        availableCPU = multiprocessing.cpu_count() 
        if availableCPU > 1:
            availableCPU -= 1
        args.threads = availableCPU
        print str(args.threads)+ ", calculated from (1 - total CPUs), on your machine will be used. To change the this, use the -n flag."
        print "More threads will run faster, but you shouldn't use more than the number of CPUs in your machine.\n"
    if args.bootstraps is None:
        print ""
        print "100 bootstraps will be performed. To change the number of bootstraps, use the -bootstraps flag.\n"
        args.bootstraps = 100
    if args.alphabet is None:
        print "No alphabet was selected, so amino acids will be used by default. Use the -alphabet flag to specify an alphabet.\n"
        args.alphabet = "AA"
    command = "python main.py " + str(args.infile) + " " + str(args.alphabet) + " " + str(args.bootstraps) + " " + str(args.threads)
    subprocess.call(command, shell=True)
    

main()
