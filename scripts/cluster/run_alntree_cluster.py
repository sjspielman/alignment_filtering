######### LAST EDIT ON 3/29/14. ##############
## Call as python run_alntree_cluster.py <run_number> <gene> <seqdir> <numproc>
## seqdir  = where the raw sim sequence files are. given in the qsub.

import re, os, sys, subprocess, shutil

sys.path.append("src/")
from aligner import *
from misc import *
from treebuilder import *
from weight_treebuilder import *
from scorer import *
from bootstrapper import *


##############################################################################################################################
############################################ INPUT ARGUMENTS FROM QSUB SCRIPT ################################################

blah = int(sys.argv[1]) - 1 ##array job so -1
gene = sys.argv[2] ## either or5, rho, prk, flat
seqdir = sys.argv[3] ## where the raw simulated files are
numproc = int(sys.argv[4]) ## number of threads

##############################################################################################################################
##############################################################################################################################


##############################################################################################################################
################################################ GUIDANCE THINGS #############################################################

n=100
refaln_file='refaln.fasta' # Will contain the reference (unmasked!) alignment
temp_res='temp.aln'        # Used as a temporary alignment file during masking	
BootDir='BootDir/'        
prepareDir(BootDir)

masks={'_30':0.3, '_50':0.5, '_70':0.7, '_90':0.9}
##############################################################################################################################
##############################################################################################################################



##############################################################################################################################
################ COPY OVER RAW SIMULATION SEQUENCES AND SET UP RETURN DIRECTORIES ############################################
## copy over raw simulation sequences.
rawaa='rawsim_aa'+str(blah)+'.fasta'
rawnuc='rawsim_codon'+str(blah)+'.fasta'
prealn_file = rawaa

command='cp '+seqdir+'/'+rawaa+' .'
call=subprocess.call(command, shell=True)
assert(call == 0), "Raw aa not copied"
command='cp '+seqdir+'/'+rawnuc+' .'
call=subprocess.call(command, shell=True)
assert(call == 0), "Raw nuc not copied"

shutil.copy(rawnuc, 'BootDir/') #need to also bring this into BootDir since will Pal2Nal at the very end there.

alndir_aa='aaguided_'+gene
alndir_nuc='nucguided_'+gene
treedir='aatrees_'+gene
miscdir='bootdirs_'+gene
os.mkdir(alndir_aa)
os.mkdir(alndir_nuc)
os.mkdir(treedir)
os.mkdir(miscdir)

finalTree_file="aatree"+str(blah)+".txt"

##############################################################################################################################
##############################################################################################################################


##############################################################################################################################
################################################ GUIDANCE MODULES ############################################################

amod = MafftAligner("/home/sjs3495/bin/bin/mafft", " --auto --quiet ")
### Note that you can align with muscle and/or clustal if you feel passionate about it, but you'll have to set this up on your own. Relevant classes in src/aligner.py 

# Tree builder (build the boostrap trees)
tmod=builderFastTree("/share/apps/fasttree-2.1.3/FastTreeMP", " -fastest -nosupport -quiet ") # -nosupport **MUST** be there

# Scoring tree. Note that raxml is having some optimization issues with protcatwag.
wtmod=weightRAxML("/share/apps/RAxML-7.7.6/bin/raxmlHPC-MPI", " -m PROTGAMMAWAG -T "+str(numproc)) # You can provide other options here if you are comfortable with RAxML.

# Scorer
smod = Scorer()

# Bootstrapper. Most things are going to happen using this class.
bmod = AllBootstrapper(bootstraps = n, prealn_file = prealn_file, refaln_file = refaln_file, BootDir = BootDir, 
                       threads = numproc, aligner=amod, tree_builder = tmod, weight_tree_builder = wtmod, scorer = smod)
##############################################################################################################################
##############################################################################################################################



##############################################################################################################################
############################################## RUNNING "GUIDANCE" HERE #######################################################

# Build reference alignment
amod.makeAlignment(prealn_file, refaln_file)

# Bootstrap. Creates perturbed guide trees and alignments and then scores according to our 6 algorithms.
(numseq, alnlen, alg_scores)=bmod.runBootstrap()	
	
# Mask across specified thresholds for all algorithms.
# Save the final masked alignment file. Can also convert to nucleotide alignment using Pal2Nal if so desired (see src/misc.py for function details)
print "\nMasking residues"
for x in masks:
	for alg in alg_scores:
		outfile = alg+x+"_"+str(blah)+".fasta"
	
		maskResiduesNOMAP(refaln_file, numseq, alnlen, alg_scores[alg], masks[x], 'fasta', temp_res, "protein")
		Pal2Nal(temp_res, rawnuc, 'fasta', 'fasta', outfile, 'fasta')
		shutil.copy(temp_res, '../'+alndir_aa+'/'+outfile)
		shutil.copy(outfile, '../'+alndir_nuc)
 
	
#Also copy reference alignment to alndir_aa and alndir+nuc.
outref='refaln'+str(blah)+'.fasta'
shutil.copy('refaln.fasta', '../'+alndir_aa+'/'+outref)
Pal2Nal('refaln.fasta', rawnuc, 'fasta', 'fasta', '../'+alndir_nuc+'/'+outref, 'fasta')

# Save the tree
command='mv scoringtree.tre ../'+treedir+'/aatree'+str(blah)+'.txt'
subprocess.call(command, shell=True)

# Save the BootDir
os.chdir('../')
new = 'bootdir'+str(blah)+'.tgz'
prepareDir(BootDir, save=True, newname=new)
shutil.copy(new, miscdir)

