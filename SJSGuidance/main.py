#!/usr/bin/python

######### GUIDANCE REIMPLEMENTATION (with some novel scoring algorithms) WRITTEN BY SJS, ETD ##########
from numpy import *
from dendropy import *

import os
import fnmatch
import re
import subprocess
import sys
import shutil
import argparse

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,generic_protein
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

sys.path.append("src/")
from aligner import *
from treebuilder import *
from weight_treebuilder import *
from scorer import *
from bootstrapper import *
from misc import *



###################### User input (or derived from user input) ###########################

n =    10      #bootstraps
numproc =  1  #threads

prefix =  "TESTSEQ" ## please set this up to be the infile's name without the extension. as in, if they have a file called "myseqs.fasta" the prefix is "myseqs"
unaligned =  "TESTSEQ.fasta" #infile
format    = "fasta" #infile format. 
alphabet  = "protein"# This should be either "protein" or "nucleotide"
final_aln_dir= "final_aln/"# where alignments will end up. name this like the infile. Make sure this string ends with /
final_boot_name = "final_boot" #this is the final name for the .tgz (compress everything in BootDir/). Can we make the name a time-stamp or something?

############################### Internal variables #######################################
prealn_file='prealn.fasta' # Will contain the raw (unaligned) sequences in fasta format with integer sequence names
refaln_file='refaln.fasta' # Will contain the reference (unmasked!) alignment
temp_res='tempaln_res.aln' # Used as a temporary alignment file during masking	

scoreTree_file='scoringtree.tre'     # Will contain the phylogeny used for computing all weights
weightfile='treeweights.txt'         # Will contain phylogenetic weights as calculated by BranchManager
dist_matrix_file = 'dist_matrix.txt' # Will contain the patristic distance matri
BootDir='BootDir/'                   # Directory where most stuff will happen

# By default, masking will happen at these four thresholds. Again, change if you want.
masks={'_30':0.3, '_50':0.5, '_70':0.7, '_90':0.9}

prepareDir(final_aln_dir)
prepareDir(BootDir)


################################# Prepare classes ########################################

### amod, tmod, wtmod all take arguments in the form (executable, options). We recommend these options, but you are welcome to play around. 

# Aligner
amod = MafftAligner("mafft", " --auto --quiet ")
### Note that you can align with muscle and/or clustal if you feel passionate about it, but you'll have to set this up on your own. Relevant classes in src/aligner.py 

# Tree builder (build the boostrap trees)
tmod=builderFastTree("FastTree", " -fastest -nosupport -quiet ") # -nosupport **MUST** be there

# But the scoring tree
if alphabet == "protein":
	model = "-m PROTCATWAG"
elif alphabet == "dna":
	model = "-m GTRCAT"
wtmod=weightRAxML("raxmlHPC", model) # You can provide other options here if you are comfortable with RAxML.

# Scorer
smod = Scorer()

# Bootstrapper. Most things are going to happen using this class.
bmod = AllBootstrapper(amod, tmod, wtmod, smod)


############################################### ACTUALLY RUN GUIDANCE HERE #################################################

# Create map for sequences
print "Initializing"
map=buildMap(unaligned, 'fasta', prealn_file)

# Build reference alignment
print "Building reference alignment"
amod.makeAlignment(prealn_file, refaln_file)

# Bootstrap. Creates perturbed guide trees and alignments and then scores according to our 6 algorithms.
(numseq, alnlen, alg_scores)=bmod.runBootstrap(BootDir, unaligned, refaln_file, n, numproc,scoreTree_file, weightfile, dist_matrix_file)	
	

# Mask across specified thresholds for all algorithms.
# Save the final masked alignment file. Can also convert to nucleotide alignment using Pal2Nal if so desired (see src/misc.py for function details)
print "\nMasking residues"
for x in masks:
	for alg in alg_scores:
		maskResidues(refaln_file, numseq, alnlen, map, alg_scores[alg], masks[x], 'fasta', temp_res, alphabet)
		outfile_aa=prefix+"_"+alg+x+".fasta"
		shutil.copy(temp_res, '../'+final_aln_dir+outfile_aa)

# Save unmasked alignment as well
unMap(map, refaln_file, '../'+final_aln_dir+prefix+"_unmasked.fasta", numseq)


# Clean up.
os.chdir('../')
os.remove(refaln_file)
os.remove(prealn_file)
prepareDir(BootDir, save=True, newname=final_boot_name)
	
	
	
	