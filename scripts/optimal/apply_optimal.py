### 6/14/13
### Apply the optimal filter to reference alignment. Residue and column-wise.
### Take the reference alignments and use the guidance algorithm to compare to the TRUE ALIGNMENT. Mask all residues with a score of 0. (keep gaps!!)
### Then build phylogenies for these.
### The scorer and masker functions have been moved into here for convenience. NEED TO COPY OVER THE GUIDANCE_SCORE PROGRAM HOWEVER!!

import re, sys, os, subprocess, shutil
from numpy import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,generic_protein
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

####################################################################################
####################################################################################

def scorer(parsed_ref, numseq, alnlen, refaa, trueaa, scorefile):

	scoreCommand='./guidance_score ' + refaa + " " + trueaa + " " + scorefile
	subprocess.call(scoreCommand, shell=True)	
	
	# List with how many gaps in each column
	numGaps = []
	for i in range(0,alnlen): #number of sequences			
		findgaps = str(parsed_ref[:,i])
		num = findgaps.count('-')
		numGaps.append(num)
	
	#Read in and normalize scores
	scores = loadtxt(scorefile)	
	normedscores=zeros(alnlen)
	final_scores = zeros(alnlen)
	t=0 #taxon counter
	s=0 #entry in row counter
	for taxon in scores:
		for score in taxon:
			norm = (numseq - numGaps[s] - 1)
			if norm==0 or score==0:
				normedscores[s]=0
			else:
				normedscores[s]=(score/norm)
			s+=1
		final_scores = vstack((final_scores, normedscores))
		normedscores=zeros(alnlen)
		s=0	
		t+=1
	final_scores = delete(final_scores, (0), axis=0)
	return final_scores

####################################################################################

def maskResidues(parsed_ref, numseq, alnlen, scores, formatout, final_file, seqType):
	
	#All things are masked to become ?
	new='?'	

	maskedMSA=MultipleSeqAlignment([])
	
	counter=0
	for row in range(numseq):
		newseq=''
		for position in range(alnlen):
			thispos=str(parsed_ref[row].seq[position])
			isgap=re.search('-', thispos)
			if isgap:
				newseq=newseq+parsed_ref[row].seq[position]
			else:
				if scores[row][position]!=1: #mask if shitty
					newseq=newseq+new
				elif scores[row][position]==1: #or, keep that position
					newseq=newseq+parsed_ref[row].seq[position]
				else:
					print "score is neither 1 nor 0, res fxn"
					break
		if str(seqType)=='protein':
			aln_record=SeqRecord(Seq(newseq,generic_protein),id='t'+str(counter), description='')
		elif str(seqType)=='nucleotide':
			aln_record=SeqRecord(Seq(newseq,generic_dna),id='t'+str(counter), description='')
		maskedMSA.append(aln_record)
		counter+=1

	outhandle=open(final_file, 'w')
	outhandle.write(maskedMSA.format(str(formatout)))
	outhandle.close()
	return 0	

####################################################################################

def maskColumns(parsed_ref, numseq, alnlen, scores, formatout, final_file, seqType):

	#All columns masked to become ?
	mcol=numseq*'?'
	raw=[]
	
	## Loop over columns in the scores.		
	col_counter=0
	for pos_score in scores.T: #transposes it so loops over columns of all_scores instead of default rows
		#print alnlen, col_counter
		# Making a bunch of lists containing indices for if it's not a gap, 'G' if it is a gap. 
		notGaps=[]
		for i in range(0,alnlen): #number of sequences			
			findgaps = str(parsed_ref[:,i])
			shortlist=[]
			for counter in range(len(findgaps)): #column!
				if findgaps[counter]!='-':
					shortlist.append(counter)
				else:
					shortlist.append('G')
			notGaps.append(shortlist)

		count=0
		sum=0
		for i in range(numseq):
			if notGaps[i] != 'G':
				sum+=pos_score[i]			
				count+=1
		colscore=float(sum)/float(count)				
		colseq=''	
		for row in range(numseq):
			colseq=colseq+parsed_ref[row].seq[col_counter]
							
		if colscore==1:
			raw.append(str(colseq))
		elif colscore!=1:
			raw.append(str(mcol))		
		else:
			print "col fxn problem with scores"
			break	
		col_counter+=1
	
	## Here, I turn them into lists of list. Each list is a ROW in the final MSA. hurray!
	raw_transposed=zip(*raw)
	## The alignments have been masked. Now need to unmap the names.
	maskedMSA=MultipleSeqAlignment([])
	counter=0
	for row in raw_transposed:
		newrow=''.join(row)
		if str(seqType)=='protein':
			aln_record=SeqRecord(Seq(newrow,generic_protein),id='t'+str(counter), description='')
		elif str(seqType)=='nucleotide':
			aln_record=SeqRecord(Seq(newrow,generic_dna),id='t'+str(counter), description='')
		maskedMSA.append(aln_record)
		counter+=1
	
	outhandle=open(final_file, 'w')
	outhandle.write(maskedMSA.format(str(formatout)))
	outhandle.close()
	return 0
		
####################################################################################
####################################################################################



## Keep first blank for indexing with array jobs.
dir=['', 'bigindel_16rc', 'medindel_16rc', 'shortindel_16rc', 'bigindel_64rc', 'medindel_64rc', 'shortindel_64rc', 'bigindel_16p', 'medindel_16p', 'shortindel_16p', 'bigindel_64p', 'medindel_64p', 'shortindel_64p', 'bigindel_16r', 'medindel_16r', 'shortindel_16r', 'bigindel_64r', 'medindel_64r', 'shortindel_64r', 'bigindel_16b', 'medindel_16b', 'shortindel_16b', 'bigindel_64b', 'medindel_64b', 'shortindel_64b']
direc=dir[int(sys.argv[1])]
aligner=str(sys.argv[2])

## Copy over all the reference files to mask and the truealn to compare to
seqdir_aa='refaln_aa'
os.mkdir(seqdir_aa)
seqdir_nuc='refaln_nuc'
os.mkdir(seqdir_nuc)
command1='cp /home/sjs3495/'+aligner+'_alntree/aaguided_'+str(direc)+'/refaln* '+seqdir_aa
call1=subprocess.call(command1, shell=True)
command2='cp /home/sjs3495/'+aligner+'_alntree/nucguided_'+str(direc)+'/refaln* '+seqdir_nuc
call2=subprocess.call(command2, shell=True)
command3='cp /home/sjs3495/'+aligner+'_alntree/aaguided_'+str(direc)+'/truealn* '+seqdir_aa
call3=subprocess.call(command3, shell=True)
command4='cp /home/sjs3495/'+aligner+'_alntree/nucguided_'+str(direc)+'/truealn* '+seqdir_nuc
call4=subprocess.call(command4, shell=True)


## Create the output directories
alndir_aa='optimal_aaguided_'+direc
alndir_nuc='optimal_nucguided_'+direc
treedir='optimal_aatrees_'+direc
os.mkdir(alndir_aa)
os.mkdir(alndir_nuc)
os.mkdir(treedir)

scorefile="scores.txt"

## 100 reference alignments to optimally mask
for n in range(100):

		# True simulated alignment. Then, the raw reference alignment from the guidance runs. THIS is what we will apply the optimal filter to.
		trueaa=seqdir_aa+'/truealn_aa'+str(n)+'.fasta'
		truenuc=seqdir_nuc+'/truealn_codon'+str(n)+'.fasta'
		refaa=seqdir_aa+'/refaln'+str(n)+'.fasta'
		refnuc=seqdir_nuc+'/refaln'+str(n)+'.fasta'
		
		#Final output files
		final_aa_res=alndir_aa+"/optres"+str(n)+".fasta"
		final_nuc_res=alndir_nuc+"/optres"+str(n)+".fasta"
		final_aa_col=alndir_aa+"/optcol"+str(n)+".fasta"
		final_nuc_col=alndir_nuc+"/optcol"+str(n)+".fasta"
		
		# Gettin some info
		infile=open(refaa, 'r')
		parsed_refaa = AlignIO.read(infile, 'fasta')
		infile.close()
		alnlen_aa=parsed_refaa.get_alignment_length()
		numseq=len(parsed_refaa)		
		infile=open(refnuc, 'r')
		parsed_refnuc = AlignIO.read(infile, 'fasta')
		infile.close()
		alnlen_nuc=parsed_refnuc.get_alignment_length()
		numseq=len(parsed_refnuc)			
		
		## Optimal filter, amino acid
		print "Scoring ",n
		scores=scorer(parsed_refaa, numseq, alnlen_aa, refaa, trueaa, scorefile)	
		maskResidues(parsed_refaa, numseq, alnlen_aa, scores, "fasta", final_aa_res, "protein")
		maskColumns(parsed_refaa, numseq, alnlen_aa, scores, "fasta", final_aa_col, "protein")

		## Optimal filter, nucleotide
		scores=scorer(parsed_refnuc, numseq, alnlen_nuc, refnuc, truenuc, scorefile)	
		maskResidues(parsed_refnuc, numseq, alnlen_nuc, scores, "fasta", final_nuc_res, "nucleotide")
		maskColumns(parsed_refnuc, numseq, alnlen_nuc, scores, "fasta", final_nuc_col, "nucleotide")
		
		
##### Build the trees.

prefix=['optres', 'optcol']
for n in range(100):
	for case in prefix:
		alnfile=alndir_aa+'/'+case+str(n)+'.fasta'
		newfile=treedir+'/'+case+str(n)+'.tre'
		command2='/share/apps/fasttree-2.1.3/FastTreeMP -nosupport -gamma -mlacc 2 -slownni -spr 4 '+alnfile+' > '+newfile
		runit=subprocess.call(command2, shell=True)