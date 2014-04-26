## Sweep across posterior probability cutoffs. Mostly useful for ROC curve.
## Note that this must be run on a gene-by-gene basis to prevent file sizes from getting out of hand.

import re, os, sys, subprocess, fnmatch, csv, shutil
from numpy import *
from Bio import AlignIO, SeqIO
import re, os, sys, subprocess, fnmatch, csv, shutil
from numpy import *
from Bio import AlignIO, SeqIO
from parsing_fxns import *
 
if len(sys.argv) != 3:
	print "Usage: python sweepfubar.py <dataset> <gene>   , where dataset is HA or GP41, gene is or5, rho, prk, or flat."
	sys.exit(0) 
 
 
datadir='/Users/sjspielman/Dropbox/aln/results/'

################### Input arguments ###################
dataset = sys.argv[1]
gene = sys.argv[2]
assert (dataset == 'HA' or dataset == 'GP41'), "Must specify either HA or GP41 as the dataset, as first arg."
assert (gene == 'or5' or gene == 'rho' or gene == 'prk' or gene == 'flat'), "Must specify gene as second arg."
if dataset == 'GP41':
	datadir += 'GP41/'
	posStart = 10
elif dataset == 'HA':
	datadir += 'HA/'
	posStart = 18
######################################################	
	

################## Important stuff ####################
# Directories: fubar output, alignments (all made with linsi, except for true alignments as generated by Indelible)
fudir   = datadir+'fubar/fubar_'+gene+'/'
alndir  = datadir+'alntree/aaguided_'+gene+'/'
	
# Directories: true simulated alignments and evolutionary rate categories
truerates_dir=datadir+'Simulation/truerates/'+gene+'/'
truealn_dir=datadir+'Simulation/sequences/'+gene+'/'

algs=['refaln', 'Guidance', 'GuidanceP', 'BMweights', 'BMweightsP', 'PDweights', 'PDweightsP']
cutoffs=arange(0,1.01,0.01)
######################################################



outfile='/Users/sjspielman/Research/alignment_filtering/data/parsed_data/revision/fubar_'+dataset+'_'+gene+'_sweep.txt'
outhandle=open(outfile, 'w')
outhandle.write('count\tcutoff\ttprate\tfprate\ttnrate\tfnrate\taccuracy\tcase\tpenal\tgene\tmethod\n')

		
for n in range(100):
	print str(n)
	
	## File names (refaln, truealn, truerates)
	refaln=alndir+'refaln'+str(n)+'.fasta'
	trfile=truerates_dir+'truerates'+str(n)+'.txt'
	truealn=truealn_dir+'truealn_aa'+str(n)+'.fasta'
			
	## Read in the reference alignment and collect some relevant info
	handle = open(refaln, 'r')
	refparsed=AlignIO.read(refaln, 'fasta')
	handle.close()
	alnlen=len(refparsed[0])
	numseq=len(refparsed)
	
	## Read in the true alignment
	handle = open(truealn, 'r')
	trueparsed=AlignIO.read(handle, 'fasta')
	handle.close()
		
	## Build map to true alignment and obtain simulated positive selection state (binary - 0=notpos, 1=pos)
	# wantRef = sites we want from reference. wantTrue = sites we want from true. Note only singletaxonmap is used.	
	wantRef, wantTrue = singleTaxonMap(trueparsed, refparsed, numseq, alnlen)	
	truepos = parseTrueRates(trfile, wantTrue, posStart)
	
	################################ Accuracy assessment ###################################
	for alg in algs:
		
		## Get file names and whether or not gap-penalized algorithm
		if alg=='refaln':
			fubar=fudir+'fubar/refaln'+str(n)+'.fasta.fubar'
			penal = 'zero'
		else:
			name = alg+'_50_'+str(n)+'.fasta'
			fubar=fudir+'fubar/'+name+'.fubar'
			if 'P' in alg:
				penal = 'yes'
			else:
				penal = 'no'
		
		# Get case info
		testprobs = parseFUBAR(wantRef, fubar)	
		assert(len(truepos)==len(testprobs)), "FUBAR Mapping has failed."
	

		## Accuracy across posterior probability cutoffs
		for x in cutoffs:
			(tp,tn,fp,fn,tprate,fprate,tnrate,fnrate,accuracy)=getAccuracy(float(x), truepos, testprobs)
			outhandle.write(str(n)+'\t'+str(x)+'\t'+str(tprate)+'\t'+str(fprate)+'\t'+str(fnrate)+'\t'+str(fnrate)+'\t'+str(accuracy)+'\t'+alg+'\t'+penal+'\t'+gene+'\tfubar\n')	
outhandle.close()

		
		
		
	
			
		
		
		
		
			
			