import re, sys, csv
from numpy import *
from Bio import AlignIO, SeqIO
import parsing

# Weighted algorithms
walgs=['BMweights', 'PDweights', 'BMweightsP', 'PDweightsP']
# Guidance algorithms
galgs=['Guidance', 'GuidanceP']
masks={'30': 'thirty', '50':'fifty', '70':'seventy', '90':'ninety'}
genes=['or5', 'rho', 'prk', 'flat']
pp_cutoff = 0.895 # Posterior probability threshold for calling sites as positively selected or not.

datadir='/Users/sjspielman/Dropbox/aln/results/'
dataset = sys.argv[1]
assert (dataset == 'HA' or dataset == 'GP41'), "Must specify either HA or GP41 as the dataset."
if dataset == 'GP41':
	datadir += 'GP41/'
	posStart = 10
elif dataset == 'HA':
	datadir += 'HA/'
	posStart = 18

outfile='/Users/sjspielman/Research/alignment_filtering/data/parsed_data/fubar_'+dataset+'_90.txt'
outhandle=open(outfile, 'w')
outhandle.write('count\ttprate\tfprate\t\tfnrate\taccuracy\tcase\tgene\tmask\tmethod\tpenal\n')


for gene in genes:

	print gene+'\n'
	
	############ Set up gene-specific data directories ############
	
	# Directories: fubar output, paml output, alignments (all made with linsi, except for true alignments as generated by Indelible)
	fudir   = datadir+'fubar/fubar_'+gene+'/'
	alndir  = datadir+'alntree/nucguided_'+gene+'/'
	
	# Directories: true simulated alignments and evolutionary rate categories
	truerates_dir=datadir+'Simulation/truerates/'+gene+'/'
	truealn_dir=datadir+'Simulation/sequences/'+gene+'/'
	
			
	for n in range(100):
		print str(n)
		
		## File names (refaln, truealn, truerates)
		refaln=alndir+'refaln'+str(n)+'.fasta'
		trfile=truerates_dir+'truerates'+str(n)+'.txt'
		truealn=truealn_dir+'truealn_codon'+str(n)+'.fasta'
		
		
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
		true_alnlen = len(trueparsed[0])
		
		###########################################################################################################
		#################### Assess accuracy for the true alignment before anything else ##########################
	
		fubar = fudir+'truealn'+str(n)+'.fasta.fubar'	
		
		# Map for truealn only. Can just go position by position as the true alignment is, shockingly, the same as itself.
		map = []
		for i in range(truealn_len):
			map.append(i)
		truepos = parseTrueRates(trfile, map, posStart)
		
		testprobs = parseFUBAR(map, fubar)	
		assert( len(truepos)==len(testprobs)), "True FUBAR Mapping has failed."
		(tp,tn,fp,fn,tprate,fprate,tnrate,fnrate,accuracy) = getAccuracy(ppcutoff, truepos, testprobs)
		outhandle.write(str(n)+'\t'+str(tprate)+'\t'+str(fprate)+'\t'+str(fnrate)+'\t'+str(accuracy)+'\ttruealn\t'+gene+'\ttrue\tfubar\ttrue\n')	
		###########################################################################################################

		###########################################################################################################
		################## Assess accuracy for the refaln. Comes first since it isn't masked.  ####################
		
		
		# Note that these values will be used for all subsequent alignments in this n rep		
		mapRef, mapTrue = consensusMap(trueparsed, refparsed, numseq, alnlen)	
		truepos = parseTrueRates(trfile, mapTrue, posStart)
		
		fubar = fudir+'refaln'+str(n)+'.fasta.fubar'	
		testprobs = parseFUBAR(mapRef, fubar)	
		assert( len(truepos)==len(testprobs)), "Reference FUBAR Mapping has failed."
		
		(tp,tn,fp,fn,tprate,fprate,tnrate,fnrate,accuracy) = getAccuracy(ppcutoff, truepos, testprobs)
		outhandle.write(str(n)+'\t'+str(tprate)+'\t'+str(fprate)+'\t'+str(fnrate)+'\t'+str(accuracy)+'\trefaln\t'+gene+'\tzero\tfubar\tzero\n')	
		###########################################################################################################

		
		###########################################################################################################		
		########################## Assess accuracy for Guidance(P), which use all masks ###########################
		for mask in masks:
			for alg in galgs:
				
				# Penalization algorithm or not? (for printing to outfile)
				if alg=='Guidance':
					penal='no'				
				else:
					penal='yes'
					
				# Collect alignment and fubar files for this algorithm
				name = alg+'_'+mask+'_'+str(n)+'.fasta'
				aln=alndir+name	
				
				# Get information relevant to this case
				parsed=AlignIO.read(aln, 'fasta')	
				fubar=fudir+name+'.fubar' 
				testprobs = parseFUBAR(mapRef, fubar)	
				assert( len(truepos)==len(testprobs)), "FUBAR Mapping has failed."
	
				## FUBAR assessment	at single posterior probability cutoff			
				(tp,tn,fp,fn,tprate,fprate,tnrate,fnrate,accuracy)=sweepRates(0.895, truepos, testprobs)
				outhandle.write(str(n)+'\t'+str(tprate)+'\t'+str(fprate)+'\t'+str(fnrate)+'\t'+str(accuracy)+'\t'+alg+'\t'+gene+'\t'+masks[mask]+'\tfubar\t'+penal+'\n')
				
		###########################################################################################################		
		####################### Assess accuracy for BM/PDweights(P), which use only mask 0.5 ######################
		for alg in walgs:		
			# Penalization algorithm or not? (for printing to outfile)
			if alg=='BMweights' or alg=='PDweights':
				penal='no'				
			else:
				penal='yes'
				
			# Collect alignment and fubar files for this algorithm
			name = alg+'_'+mask+'_'+str(n)+'.fasta'
			aln=alndir+name	
			
			# Get information relevant to this case
			parsed=AlignIO.read(aln, 'fasta')	
			fubar=fudir+name+'.fubar' 
			testprobs = parseFUBAR(mapRef, fubar)	
			assert( len(truepos)==len(testprobs)), "FUBAR Mapping has failed."
	
			## FUBAR assessment	at single posterior probability cutoff			
			(tp,tn,fp,fn,tprate,fprate,tnrate,fnrate,accuracy)=sweepRates(0.895, truepos, testprobs)
			outhandle.write(str(n)+'\t'+str(tprate)+'\t'+str(fprate)+'\t'+str(fnrate)+'\t'+str(accuracy)+'\t'+alg+'\t'+gene+'\tfifty\tfubar\t'+penal+'\n')		


outhandle.close()

			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
		
				
			
			
			
			
			
			