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


n = 100  # bootstrap n times
prealn_file='prealn.fasta'
refaln_file='refaln.fasta'
weightfile='treeweights.txt'
dist_matrix_file = 'dist_matrix.txt'
scoreTree_file='scoringtree.tre'
BootDir='BootDir/'
numproc = 2


#Final output files
simcount=1
final_guidance="guidance"+str(simcount)+".fasta"
final_BMweights="BMweights"+str(simcount)+".fasta"
final_PDweights="PDweights"+str(simcount)+".fasta"

final_guidance_p="guidance_p"+str(simcount)+".fasta"
final_BMweights_p="BMweights_p"+str(simcount)+".fasta"
final_PDweights_p="PDweights_p"+str(simcount)+".fasta"


finalscore_fileG="scoresG"+str(simcount)+".txt"
finalscore_fileBM="scoresBM"+str(simcount)+".txt"
finalscore_filePD="scoresPD"+str(simcount)+".txt"

finalscore_fileG_p="scoresG_p"+str(simcount)+".txt"
finalscore_fileBM_p="scoresBM_p"+str(simcount)+".txt"
finalscore_filePD_p="scoresPD_p"+str(simcount)+".txt"

final_BStrees_file="BStrees"+str(simcount)+".txt"
finalTree_file="aatree"+str(simcount)+".txt"


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

########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################

# Create map for sequences (required for mafft aligner, but can keep for all. doesn't waste time.)
map=mapmod.ids2int(unaligned, 'fasta', prealn_file)
	
# Build reference alignment
amod.makeAlignment(prealn_file, refaln_file)
shutil.copy(refaln_file, BootDir)

# Bootstrap
(numseq, alnlen, gscores, bmscores, pdscores, gscores_p, bmscores_p, pdscores_p)=bmod.runBootstrap(BootDir, raw, refaln_file, n, numproc, finalscore_fileG, finalscore_fileBM, finalscore_filePD, finalscore_fileG_p, finalscore_fileBM_p, finalscore_filePD_p, scoreTree_file, weightfile, dist_matrix_file)	
	
#residue masking and pal2nal

masks={'30_':float(0.3)}
algs={'guidance':guidance, 'BMweights':BMweights, 'PDweights':PDweights, 'guidance_p':guidancep, 'BMweights_p':BMweightsp, 'PDweights_p':PDweightsp}

for x in masks:
	for alg in algs:
		outfile=alg+x+str(simcount)+".fasta"
		mmod.maskResidues(refaln_file, numseq, alnlen, algs[alg], masks[x], 'fasta', temp_res, "protein", simcount, save_x_file, alg)
		Pal2Nal(temp_res, rawnuc_ints, 'fasta', 'fasta', outfile, 'fasta')
		shutil.copy(outfile, '../'+alndir_nuc)
		shutil.copy(temp_res, '../'+alndir_aa+'/'+outfile)
	
	
#Also copy reference alignment to alndir_aa and alndir+nuc.
outref='refaln'+str(simcount)+'.fasta'
shutil.copy('refaln.fasta', '../'+alndir_aa+'/'+outref)
Pal2Nal('refaln.fasta', rawnuc_ints, 'fasta', 'fasta', '../'+alndir_nuc+'/'+outref, 'fasta')






# Clean up BootDir. compress????
os.chdir('../')
bootfiles=os.listdir(BootDir)
for file in bootfiles:
	os.remove(BootDir+file)	