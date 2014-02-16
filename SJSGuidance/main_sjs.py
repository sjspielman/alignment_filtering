#### Python script to run the guidances.
#### Input files must be in fasta format. Scroll down to see settings...

#!/usr/bin/python
import sys
import os
import fnmatch
import re
import subprocess
import shutil
import argparse

sys.path.append("src/")

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


#####################################################################################
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

def prepareDir(directory, save=False, newname=None):
	''' Double check that BootDir exists and clear it of all files, as needed. ''' 
	if save:
		# Be sure newname is properly named. Yes, not an ideal way. I know. But it's fine.
		if ".tgz" not in newname:
			newname+=".tgz"
		compress = "tar -czf "+newname+" "+directory
		subprocess.call(compress, shell=True)
	
	# Cleanup and/or make
	if (os.path.exists(directory)):
		bootfiles=os.listdir(directory)
		for file in bootfiles:
			os.remove(directory+file)	
	else:
		os.mkdir(directory)
		

######################################################################################
######################################################################################


n = 1 # bootstrap n times
unaligned='TESTSEQ.fasta'
prealn_file='prealn.fasta'
refaln_file='refaln.fasta'
weightfile='treeweights.txt'
dist_matrix_file = 'dist_matrix.txt'
scoreTree_file='scoringtree.tre'

numproc = 2

final_aln_dir='alignments/'
BootDir='BootDir/' ## KEEP

prepareDir(final_aln_dir)
prepareDir(BootDir)

#Final output files. Alignment file names defined later as algorithm_maskthreshold_simcount_nuc/aa.fasta (eg, guidance50_1_nuc.fasta or guidance_p50_1_nuc.fasta)
simcount=1
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
wtmod=weightRAxML("raxmlHPC", " -m PROTCATWAG ")	

mapmod = Map()
smod = Scorer()
bmod = AllBootstrapper(amod, tmod, wtmod, smod)
mmod = Masker(bmod)

########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################


# Create map for sequences (required for mafft aligner, but can keep for all. doesn't waste time.)
map=mapmod.ids2int(unaligned, 'fasta', prealn_file)
	
# Build reference alignment
amod.makeAlignment(prealn_file, refaln_file)

# Bootstrap
(numseq, alnlen, gscores, bmscores, pdscores, gscores_p, bmscores_p, pdscores_p)=bmod.runBootstrap(BootDir, unaligned, refaln_file, n, numproc, finalscore_fileG, finalscore_fileBM, finalscore_filePD, finalscore_fileG_p, finalscore_fileBM_p, finalscore_filePD_p, scoreTree_file, weightfile, dist_matrix_file)	
	
# Residue masking and pal2nal. Can mask at multiple cutoffs as desired here.
masks={'30_':float(0.3)}
algs={'guidance':gscores, 'BMweights':bmscores, 'PDweights':pdscores, 'guidance_p':gscores_p, 'BMweights_p':bmscores_p, 'PDweights_p':pdscores_p}
temp_res='tempaln_res.aln'	

print "\nMasking residues"
for x in masks:
	for alg in algs:
		mmod.maskResidues(refaln_file, numseq, alnlen, algs[alg], masks[x], 'fasta', temp_res, "protein")
		
		# Save the final masked alignment file. Can also convert to nucleotide alignment using Pal2Nal if so desired.
		outfile_aa=alg+x+str(simcount)+"_aa.fasta"
		shutil.copy(temp_res, '../'+final_aln_dir+outfile_aa)
		#outfile_nuc=alg+str(simcount)+"_nuc.fasta"
		#Pal2Nal(temp_res, rawnuc_ints, 'fasta', 'fasta', outfile_nuc, 'fasta')
		
	
# Save unmasked alignment as well
outref='refaln'+str(simcount)+'.fasta'
shutil.copy('refaln.fasta', '../'+final_aln_dir)
#Can convert this to nucleotide if desired with Pal2Nal.


# Clean up BootDir.
os.chdir('../')
prepareDir(BootDir, save=True, newname="savedirectory") #save=True will compress and save the directory to whatever we make newname
	
	
	
	
	