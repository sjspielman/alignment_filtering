### 8/5/13. FOR USE ON CLUSTER. NO MAPPING!!!! sequences were simulated to have int leaf names.
# Currently, FastTree for bootstrap and RAxML for weighting tree. Performs unweighted (guidance), branch manager weighted (BMweights), patristic-weighted (PDweights).
# Aligner can be toggled based on sys.argv 5

import re, os, sys, subprocess, shutil
from aligner import *
from treebuilder import *
from weight_treebuilder import *
from scorer import *
from bootstrapper import *
from treebuilder import *
from numpy import *
from dendropy import *
from masker import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,generic_protein
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment


################################################################################################
################################################################################################
def Pal2Nal(palfile, nucfile, paltype, nuctype, outfile, outputformat):

	aln_parsed=list(SeqIO.parse(str(palfile), str(paltype)))
	nuc_parsed=list(SeqIO.parse(open(nucfile, 'rU'), str(nuctype)))
	
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
				nal=nal+'---'
			elif position=='?':
				nal=nal+'???'
				start+=3
				end+=3
			elif position=='X':
				nal=nal+'NNN'
				start+=3
				end+=3
			#If amino acid there, append corresponding codon
			else:
				codon=str(nuc_seq[start:end])
				nal=nal+codon
				start+=3
				end+=3
		#Make nucleotide MSA object
		aln_record=SeqRecord(Seq(str(nal), generic_dna), id=pal_id, description='')
		nucMSA.append(aln_record)
	#write alignment to file
	outfile=open(outfile, 'w')
	umm=AlignIO.write(nucMSA, outfile, str(outputformat))
	outfile.close()
	return 0
	

################################################################################################
################################################################################################


simcount = int(sys.argv[1])
if simcount == 100:  ## stupid -t
	simcount = 0
gene = sys.argv[2]
direc = sys.argv[3]
rdir = sys.argv[4]
aligner = sys.argv[5] ## mafft, clustal, or muscle


## copy over raw simulation sequences.
raw='rawsim_aa'+str(simcount)+'.fasta'
rawnuc='rawsim_codon'+str(simcount)+'.fasta'
command='cp -r /home/sjs3495/rawdata_HA_b/'+gene+'/'+direc+'/'+raw+' .'
call=subprocess.call(command, shell=True)
command='cp -r /home/sjs3495/rawdata_HA_b/'+gene+'/'+direc+'/'+rawnuc+' .'
call=subprocess.call(command, shell=True)
shutil.copy(raw, 'BootDir/')
shutil.copy(rawnuc, 'BootDir/') #need to also bring this into BootDir since will Pal2Nal at the very end there.


alndir_aa='aaguided_'+aligner+'_'+gene+'_'+direc
alndir_nuc='nucguided_'+aligner+'_'+gene+'_'+direc
treedir='aatrees_'+aligner+'_'+gene+'_'+direc
miscdir='misc_'+aligner+'_'+gene+'_'+direc
os.mkdir(alndir_aa)
os.mkdir(alndir_nuc)
os.mkdir(treedir)
os.mkdir(miscdir)



#Guidance files

n = 100  # bootstrap n times
prealn_file='prealn.fasta'
refaln_file='refaln.fasta'
weightfile='treeweights.txt'
dist_matrix_file = 'dist_matrix.txt'
scoreTree_file='scoringtree.tre'
BootDir='BootDir/'
save_x_file = rdir+'/xfile_'+aligner+'_'+gene+'_'+direc+'.txt'


#Guidance modules

if aligner == 'linsi':
	amod = MafftAligner("/home/sjs3495/bin/bin/mafft", " --auto --quiet " ) ## can add auto in for better maffting
elif aligner == 'clustal':
	amod = ClustalAligner("/home/sjs3495/bin/clustalw2/clustalw2", " -quiet ")
elif aligner == 'muscle':
	amod = MuscleAligner("/home/sjs3495/bin/muscle", " -quiet ")
else:
	print "bad aligner!!"
	assert 1==0
numproc = 4
tmod=builderFastTree("/share/apps/fasttree-2.1.3/FastTreeMP", " -quiet -nosupport ") 
wtmod=bothRAxML("/share/apps/raxmlHPC-7.3.0/bin/raxmlHPC", " -m PROTCATWAG -p 35325 ")	
smod = Scorer()
bmod = AllBootstrapper(amod, tmod, wtmod, smod)
mmod = Masker(bmod)





#Final output files
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

##############################################################################################################
##################################### RUNNING "GUIDANCE" HERE ################################################
print "Starting the Guidances"

#Build initial MSA
amod.makeAlignment(raw, refaln_file)
shutil.copy(refaln_file, BootDir)

#SCORING.
(numseq, alnlen, gscores, bmscores, pdscores, gscores_p, bmscores_p, pdscores_p)=bmod.runBootstrap(BootDir, raw, refaln_file, n, numproc, finalscore_fileG, finalscore_fileBM, finalscore_filePD, finalscore_fileG_p, finalscore_fileBM_p, finalscore_filePD_p, scoreTree_file, weightfile, dist_matrix_file)	
	
###### MASKING SCORES HERE ########################
print "Masking residues"

temp_res='tempaln_res.aln'	

###### MASKING SCORES HERE ########################
print "Masking residues"

temp_res='tempaln_res.aln'	

masks={'30_':float(0.3), '50_':float(0.5), '70_':float(0.7), '90_':float(0.9)}
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


## Save the bootstrap trees, the RAxML tree that was made during weighting, and the file of final scores
command='mv BStrees.tre ../'+miscdir+"/"+final_BStrees_file
subprocess.call(command, shell=True)
command='mv '+scoreTree_file+' ../'+treedir+'/'+finalTree_file
subprocess.call(command, shell=True)
command='mv '+finalscore_fileG+' ../'+miscdir
subprocess.call(command, shell=True)
command='mv '+finalscore_fileBM+' ../'+miscdir
subprocess.call(command, shell=True)
command='mv '+finalscore_filePD+' ../'+miscdir
subprocess.call(command, shell=True)
command='mv '+finalscore_fileG_p+' ../'+miscdir
subprocess.call(command, shell=True)
command='mv '+finalscore_fileBM_p+' ../'+miscdir
subprocess.call(command, shell=True)
command='mv '+finalscore_filePD_p+' ../'+miscdir
subprocess.call(command, shell=True)