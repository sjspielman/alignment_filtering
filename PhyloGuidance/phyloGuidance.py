import subprocess
import argparse
import sys
import os
try:
    from dendropy import *
except:
    print "Please install Dendropy. See the README for details."


def parse_args():
    parser = argparse.ArgumentParser(prefix_chars = '+-', usage = 'Enter python phyloGuidance.py --help/-h for instructions')
    parser.add_argument("-infile", dest = "infile", default = None, type = str, help = "A file containing unaligned sequences in FASTA format. CAUTION: no sanity checking performed for this!")
    parser.add_argument("-cpu", dest = "threads", default = 1, type = int, help="Number of processes to use")
    parser.add_argument("-bootstraps", dest="bootstraps", default = 100, help = "The number of bootstraps to perform.")
    parser.add_argument("-alphabet", dest = "alphabet", default = "protein", type = str, help = "Alignment alphabet. Either protein or dna.")

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
    	print
        args.infile = raw_input("Provide an input sequence file in FASTA format.\n CAUTION - no sanity checking performed for file type!\n: ")
    while not os.path.exists(str(args.infile)):
    	print
    	args.infile = raw_input("Provided input file does not appear to exist. Try again? \n")
    if args.threads is None:
        print
        import multiprocessing 
        availableCPU = multiprocessing.cpu_count() 
        if availableCPU > 1:
            availableCPU -= 1
        args.threads = availableCPU
        print str(args.threads)+ ", calculated from (1 - total CPUs), on your machine will be used. To change the this, use the -cpu flag."
        print "More threads will run faster, but you shouldn't use more than the number of CPUs in your machine.\n"
    if args.bootstraps is None:
        print
        print "100 bootstraps will be performed. To change the number of bootstraps, use the -bootstraps flag.\n"
        args.bootstraps = 100
    if args.alphabet is None:
    	print
        print "No alphabet was selected, so amino acids will be used by default. Use the -alphabet flag to specify an alphabet (dna or protein)_.\n"
        args.alphabet = "protein"
        
    print "Now running PhyloGuidance\n" 
    command = "python src/run_phyloGuidance.py " + str(args.infile) + " " + str(args.alphabet) + " " + str(args.bootstraps) + " " + str(args.threads)
    subprocess.call(command, shell=True)
    

main()
