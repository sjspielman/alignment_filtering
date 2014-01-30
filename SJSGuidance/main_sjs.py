#### Python script to run the guidances.
#### Input files must be in fasta format. Scroll down to see settings...

#!/usr/bin/python
from aligner import *
from treebuilder import *
from weight_treebuilder import *
from scorer import *
from bootstrapper import *
from map import *
from masker import *
from dendropy import *

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,generic_protein
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

from numpy import *
import os
import fnmatch
import re
import subprocess
import sys
import shutil
import argparse


######################################################################################
######################################################################################
def Pal2Nal(palfile, nucfile, paltype, nuctype, outfile, outputformat):
	''' Convert a protein alignment to a nucleotide alignment. Can handle only the ambiguities N, X, ?. '''
	## Arguments:
	##	-> palfile = protein alignment file
	##  -> nucfile = unaligned nucleotide file. Sequences should have same name and be in same order as in palfile
	##  -> paltype = palfile format (eg fasta...)
	##  -> nuctype = nucfile format
	##  -> outfile = output file for nucleotide alignment
	##  -> outputformat = format for output nucleotide alignment
	
	
	aln_parsed=list(SeqIO.parse(str(palfile), str(paltype)))
	nuc_parsed=list(SeqIO.parse(str(nucfile), str(nuctype)))	
	
	if len(aln_parsed)!=len(nuc_parsed):
		print palfile+' '+nucfile+' have different number of sequences! Please make sure that these two files correspond and that all stop codons are removed.'
		assert 1==0
	else:
		numseq=len(aln_parsed)
	
	nucMSA=MultipleSeqAlignment([])
	for p in range(0,numseq):
		pal_seq=str(aln_parsed[p].seq) #aa alignment sequence
		pal_id=str(aln_parsed[p].id)
		for n in range(0,numseq):
			if nuc_parsed[n].id==pal_id:
				nuc_seq=str(nuc_parsed[n].seq)
				nal=str()
				start=0 #counter for codon starting position
				end=3 #counter for codon ending position
				for position in pal_seq:
					#If gapped, missing, or masked position in alignment, append 3 gaps/missing/NNN to new string 
					if position=='-':
						codon='---'
						nal=nal+codon
					elif position=='?':
						codon='???'
						nal=nal+codon
						start+=3
						end+=3
					elif position=='X':
						codon='NNN'
						nal=nal+codon
						start+=3
						end+=3
					#If amino acid there, append corresponding codon
					else:
						codon=str(nuc_seq[start:end])
						nal=nal+codon
						start+=3
						end+=3
				#Make nucleotide MSA object
				nal_seq=Seq(nal)
				aln_record=SeqRecord(Seq(nal, generic_dna), id=pal_id, description='')
				nucMSA.append(aln_record)
			else:
				continue
	#write alignment to file
	outfile=open(outfile, 'w')
	umm=AlignIO.write(nucMSA, outfile, "fasta")
	outfile.close()
	return 0
######################################################################################
######################################################################################


### FILES ###
######### !!!!!! DO NOT CHANGE THESE. Or do, but no need. It won't cause a problem but really just leave these as they are. ##
BootDir='BootDir/'
prealn_file='prealn.fasta' 
refaln_file='refaln.fasta'
##########


#User files
unaligned='protein.fasta' ## Relative path to file you want to file. Should contain unaligned sequences in FASTA FORMAT.
seqType='protein' ## protein or nucleotide, but i guess everything is protein?
#rawnuc="../rawsim_codon0.fasta"
#save_x_file='../savexfile.txt'

#User options (currently set as default)
#n = int(sys.argv[1])  #number of bootstraps
#if n=='':
#	n=100
n=10
threshold=0.90 # scoring cutoff for masking residues and/or columns. keep only >=x
nproc=2 ##  1 processor default.
pflag = 0; # First, I'm hilarious. Second, 0=no gap penalization, 1=gap penalization.

finalaln_file_nuc = 'aln_nuc.fasta'
finalaln_file_aa  = 'aln_aa.fasta'
##########################################################################################
###########   DON'T CHANGE ANY OF THE FOLLOWING IN HERE. EVER. I'M SERIOUS.  #############

### each module is (executable, options) ####

# Aligner
amod = MafftAligner("mafft", " --quiet ")
###amod = MuscleAligner("muscle", " -quiet ")
###amod = ClustalAligner("clustalw2", " -quiet ")

# Tree builder (bootstrap trees)
tmod=builderFastTree("FastTree", " -fastest -nosupport -quiet ") # -nosupport MUST be there
###tmod=builderSemphy("../semphy/semphy", " -a 20 --jtt -H -J -v 5 --BPrepeats=100 ") ## MUST BE -v 5

# if weighted algorithm, wtmod needed
wtmod=scoreTreeRAxML("raxmlHPC", " -m PROTCATWAG ")	
mapmod = Map()
smod = Scorer()
bmod = Bootstrapper(amod, tmod)
mmod = Masker(bmod)

BootDir="BootDir/"

# Create map for sequences (required for mafft aligner, but can keep for all. doesn't waste time.)
map=mapmod.ids2int(unaligned, 'fasta', prealn_file)
	
# Build reference alignment
amod.makeAlignment(prealn_file, refaln_file)

# Bootstrap
(numseq, alnlen) = bmod.bootstrap(BootDir, prealn_file, refaln_file, n, numprocesses)


## something about pflag????
if Guidance:
	score_file = 'Guidance_scores.txt'
	gscores = scoreMSA_Guidance(refaln_file, n, numseq, alnlen, pflag, score_file)
	mmod.maskResidues(refaln_file, numseq, alnlen, gscores, threshold, idmap, formatout, finalaln_file, seqType)
else:
	weightTree_file = 'weighting_tree.tre'
	wtmod.buildScoreTree(refaln_file, weightTree_file)
	
	if BMweights:
		score_file = 'BMweights_scores.txt'
		bm_file = 'BranchManager_weights.txt'
		ordered_weights = wtmod.calcBMweights(weightTree_file, bm_file, numseq)
		bmscores = scoreMSA_BMweights(refaln_file, n, numseq, alnlen, pflag, ordered_weights, bm_file, score_file)
		
	if PDweights:
		score_file = 'PDweights_scores.txt'
		pd_file = 'patristic_distance_matrix.txt'
		dist_matrix = wtmod.calcPDweights(weightTree_file, pd_file, numseq)
		pdscores = scoreMSA_PDweights(refaln_file, n, numseq, alnlen, pflag, dist_matrix, pd_file, score_file)








#residue masking and pal2nal
totalmasked=mmod.maskResidues(refaln_file, numseq, alnlen, gscores, x, map, 'fasta', finalaln_file_aa, 'protein', i, save_x_file, 'pat')

#Pal2Nal(finalaln_file_aa, rawnuc, 'fasta', 'fasta', finalaln_file_nuc, 'fasta')


# Clean up BootDir. compress????
os.chdir('../')
bootfiles=os.listdir(BootDir)
for file in bootfiles:
	os.remove(BootDir+file)	

	





def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("protein_file", help="A file containing unaligned AA sequences in FASTA format", required=False, dest="infile", type=str)
    parser.add_argument("protein_format", help="Specifies the format of the input file (fasta or nex", dest=form, type=str, default="fasta")
    parser.add_argument("num_procs", type=int, help="Number of processes to use", default=1)
    parser.add_argument("bootstraps", help="The number of bootstraps to perform", required=False,
    					dest="bootstraps", default=10)
    parser.add_argument("alphabet", help="Whether AAs or NTs are used", type=str,
    					default="AA", required=False) ##AA or NT, default is AA
    parser.add_argument("gap_penalization", help="Type of gap penalization", default=0,
    					type=int, dest="pflag")
    args = parser.parse_args()
