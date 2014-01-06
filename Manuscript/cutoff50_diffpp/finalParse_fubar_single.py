## 8/1/13 SJS. Consensus mapping for assessing tprate etc.
############ ok, final redo on 12/13/13. going to do SINGLE at these cutoffs: 
## fubar. or5: 0.7083, rho:0.748, prk:0.785, flat:0.7855
## paml.  or5: 0.518, rho:0.0.3025, prk:0.0625, flat:N/A


import re, os, sys, subprocess, fnmatch, csv, shutil
from numpy import *
from Bio import AlignIO, SeqIO

################################################################################################
################################################################################################
def parseTrueRates(trfile, truealn, mapTrue):
	'''Retrieve data for each position in the true alignment, from a file generated during simulation giving the TRUE SIMULATED rates. Info starts at line 11 of the truerates files.'''
	'''ONLY FOR SITES WHERE REFTAXON IS NOT A GAP IN THE TRUEALN!!!!'''
	poslist=[] ## each entry corresponds to a position. 0=negative, 1=positively selected.
	
	## FOR HA, CAT >=18 IS POSITIVE	

	## Parse truerates file
	infile=open(trfile, 'r')
	truelines=infile.readlines()
	infile.close()
	truelines=truelines[10:] ## only keep these lines since before that it's all header crap.	
	for counter in mapTrue:
		find=re.search('^\d+\t(\d+)\t', truelines[counter])
		if find:
			rate=int(find.group(1))
			if rate >=18: ## HA!!!
				poslist.append(rate)
			else:
				poslist.append(0)
		else:
			print "bad find"
	return poslist
################################################################################################
def buildMap(trueparsed, parsed, numseq, alnlen):
	''' Returns, for each taxon in the alignment, the position that each codon has in the true alignment. Uses numpy arrays.'''
	codons=alnlen/3
	allMaps = zeros(codons, dtype=int)	## numseq and alnlen refer to parsed (NOT the true alignment!!)	
	
	taxonlist=[]
	for entry in parsed:
		taxonlist.append(str(entry.id))	
	
	for reftaxon in taxonlist:
		truelist=[]
		map=zeros(codons, dtype=int)	
		# Get sequences for this reftaxon
		for entry in trueparsed:
			if str(entry.id)==reftaxon:
				trueseq=str(entry.seq)
				break
		for entry in parsed:
			if str(entry.id)==reftaxon:
				refseq=str(entry.seq)
				break
		#Build the map. Record the index for each non-gap site in trueseq. 
		#Then, go through the refseq and for the nongaps, pop off that index. For the gaps, add a 'G'.
		for a in range(0,len(trueseq),3):
			if trueseq[a] != '-': 
				truelist.append(a/3)		
		for a in range(0,len(refseq),3):
			if refseq[a] == '-':
				map[a/3] = '-1'
			else:
				map[a/3]=truelist.pop(0)
		printmap=map.tolist()
		#print printmap,'\n'
		allMaps = vstack((allMaps, map))
	allMaps = delete(allMaps, 0, axis=0)	

	printall = allMaps.tolist()
	#print allMaps

	# Now get the consensus
	(mapRef, mapTrue) = getConsensus(allMaps)
	
	return (mapRef, mapTrue)
	
	
def getConsensus(allMaps):
	consMap=[] # will have final consensus map! hurray.
	for column in allMaps.T: # transpose so can loop over columns.
		# Check that the most frequent occurs at least 50% of the time.			
		col = column.tolist()
		most = int(max(set(col), key=col.count))
		countmost = col.count(most)
		perc_max = float(countmost)/float(len(col))
		if perc_max >= 0.5:
			consMap.append(most)
		else:
			consMap.append(-1)	
	#print consMap
	# Create two lists to tell which indices to collect from reference aln and which from true aln.
	mapRef=[]
	mapTrue=[]
	counter=0
	for entry in consMap:
		if entry != -1:
			mapRef.append(counter)
			mapTrue.append(entry)
		counter+=1	
	return (mapRef, mapTrue)
################################################################################################
def parseFubar(map, fufile):
	'''For the relevant positions, retrieve the pr(alpha>beta)=pr(positively selected). Based on INDEX in the map.'''
	## Read in fubar file
	testprobs=[] ## contains the prob(alpha>beta) values for the truefubar results
	fubar=csv.reader(open(fufile,'r'))
	allfubar=[]
	
	for row in fubar:
		if row[0]=='Codon':
			continue
		else:
			allfubar.append(float(row[4]))
	for entry in map:
		testprobs.append(allfubar[entry])	
	return testprobs
################################################################################################

def sweepRates(x, truepos, testprobs):
	'''This function doesn't actually sweep. Just compares probabilities for a given posterior probability, x.'''
	tp=0
	tn=0
	fp=0
	fn=0
	x=float(x)

	for i in range(len(truepos)):
		#Positive cases
		if float(testprobs[i])>=x:
			if truepos[i]==0:
				fp+=1
			elif truepos[i]>0:
				tp+=1
				#if str(x)=='0.9':
				#	print "tp", truepos[i]
		#Negative cases
		elif float(testprobs[i])<x:
			if truepos[i]==0:
				tn+=1
			elif truepos[i]>0:
				fn+=1
				#if str(x)=='0.9':
				#	print "fn", truepos[i]
	(tprate,fprate,tnrate,fnrate,accuracy) = calcStats(tp, fp, tn, fn)
	return (tp,tn,fp,fn,tprate,fprate,tnrate,fnrate,accuracy)

def calcStats(tp, fp, tn, fn):
	'''Pretty obvious.'''
	if fn==0 and tp==0:
		fnrate=0
		tprate=0
	else:
		fnrate=(float(fn))/(float(fn+tp))
		tprate=(float(tp))/(float(fn+tp))
	if tn==0 and fp==0:
		tnrate=0
		fprate=0
	else:
		tnrate=(float(tn))/(float(tn+fp))
		fprate=(float(fp))/(float(tn+fp))
	accuracy = (float(tp+tn)/float(tp+tn+fp+fn))
	return (tprate,fprate,tnrate,fnrate,accuracy)
################################################################################################
def percentMaskedRes(parsed, numseq, alnlen):
	'''Collects information about gappiness and masking.'''
	# For masking: (total number of ?) / (total number of letters plus ?). This way, not an average per taxon.
	taxamasked=0     # number of taxa where any residues have been masked
	totalmasked=0    # total number of gaps in the entire alignment
	totalgaps=0      # total number of gaps in the entire alignment
	numres=0
	for entry in parsed:
		seq=str(entry.seq)
		nummasked = seq.count('?')
		if nummasked!=0:
			taxamasked+=1
		totalmasked=totalmasked+nummasked
		numgaps = seq.count('-')
		totalgaps=totalgaps+numgaps
		numres=numres+len(seq) - numgaps
	#Percent masked/gaps in alignment
	percentmasked=float(totalmasked)/float(numres)				## PERCENT RESIDUES IN WHOLE ALN
	percentgaps=float(totalgaps)/(float(numseq)*float(alnlen)) ## PER SEQUENCE AVERAGE	
	return (taxamasked, percentmasked, percentgaps)
################################################################################################
def findColMasking(parsed, alnlen, numseq, map):
	'''Tells you how masked a column is, provided it is in the consensus map. However apparently this isn't being used currently.'''
	percmask=[] # a percentage for each column for how many RESIDUES are masked out of num residues in column (ignore gaps!!)
	allcounts=[]
	truth=[]
	
	seqs=[]
	for entry in parsed:
		seqs.append(str(entry.seq))
		
	for a in range(0, alnlen, 3):
		if a/3 in map:
			col=''
			for s in range(numseq):
				col=col+str(seqs[s][a]) # don't need whole codon since first position tells what everything else is anyways.
			numGaps = col.count('-')
			countmask = col.count('?')
			countres = numseq - numGaps
			perc=float(countmask)/float(countres)
			percmask.append(perc)
			allcounts.append(countmask)
	return(percmask, allcounts)		

def compareProbs(ref, test, true):
	'''Compares the different inferences at each site between the reference and the masking cases. Returns the net difference in those values!!'''
	
	## As of now, the ref and test contain probabilities rather than binary @ 0.9. Convert them
	ref2=[]
	test2=[]
	for entry in ref:
		if float(entry) >= 0.9:
			ref2.append(1)
		else:
			ref2.append(0)
	for entry in test:
		if float(entry)>=0.9:
			test2.append(1)
		else:
			test2.append(0)
	tp=0
	fp=0
	tn=0
	fn=0
	difsites=0
	for a in range(len(ref2)):
		if ref2[a] != test2[a]:
			#print ref2[a], test2[a], true[a]
			difsites+=1
			if true[a] == 1 and ref2[a] == 1:
				tp-=1
				fn+=1
			elif true[a] == 1 and ref2[a] == 0:
				tp+=1
				fn-=1	
			elif true[a] == 0 and ref2[a] == 1:
				tn+=1
				fp-=1
			elif true[a] == 0 and ref2[a] == 0:
				tn-=1
				fp+=1
	return(difsites, tp, fp, tn, fn)
################################################################################################
################################################################################################
################################################################################################
################################################################################################


prefix=['guidance', 'BMweights', 'PDweights', 'guidance_p', 'BMweights_p', 'PDweights_p']
prefs=['refaln'] # can also contain truealn 
genes={'or5':0.7083, 'rho':0.748, 'prk':0.785, 'flat':0.7855} ## these are the fubar pp cutoffs which yield a fpr=1%
masks={'50':'fifty', '70':'seventy'}
base='seqs_real'


outfile='analysis/fubarsingle.txt'
outhandle=open(outfile, 'w')
outhandle.write('count\ttprate\tfprate\ttnrate\tfnrate\taccuracy\ttaxamasked\tpercentmasked\tpercentgaps\tcase\tgene\taligner\tmethod\tmask\n')


for gene in genes:

	print gene+'\n'

	fudir = 'fubar/fubar_'+gene+'_'+base+'/'
	alndir = 'alntree/nucguided_linsi_'+gene+'_'+base+'/'
	
	truerates_dir='../Simulation/truerates/'+gene+'/'+base+'/'
	truealn_dir='../Simulation/sequences/'+gene+'/'+base+'/'
	#truefu_dir='truefu/truefu_'+gene+'_'+base+'/'
	
			
	for n in range(100):
		print str(n)
		
		trfile=truerates_dir+'truerates'+str(n)+'.txt'
		truealn=truealn_dir+'truealn_codon'+str(n)+'.fasta'
		#truefu=truefu_dir+'truealn'+str(n)+'.fubar'
		
		refaln=alndir+'refaln'+str(n)+'.fasta'
		reffu=fudir+'refaln'+str(n)+'.fasta.fubar'
		
		trueparsed=AlignIO.read(truealn, 'fasta')
		refparsed=AlignIO.read(refaln, 'fasta')
		alnlen=len(refparsed[0])
		numseq=len(refparsed)
		
		# Build map to true alignment. This map will serve for all cases associated w/ this reference alignment.
		(mapRef, mapTrue) = buildMap(trueparsed, refparsed, numseq, alnlen)
		
		# Get the true rates (binary 1/0 for pos/not pos)
		truepos =  parseTrueRates(trfile, truealn, mapTrue)

				
		for case in prefs:
		
			if case == 'refaln':
				casealn=refaln
				casefu=reffu
				parsed=refparsed
				testprobs = parseFubar(mapRef, casefu)	
				
			if case == 'truealn':
				parsed=trueparsed
				testprobs = parseFubar(mapTrue, truefu)
				
			## Ensure that the mapping went ok.
			if len(truepos)!=len(testprobs):
				print case, len(map), len(truepos), len(testprobs)
				print "There is a serious flaw in your logic. Go cry in the corner.\n\n"
				assert 1==0
				
			(tp,tn,fp,fn,tprate,fprate,tnrate,fnrate,accuracy)=sweepRates(float(genes[gene]), truepos, testprobs)
			(taxamasked, percentmasked, percentgaps)=percentMaskedRes(parsed, alnlen, numseq)
			outhandle.write(str(n)+'\t'+str(tprate)+'\t'+str(fprate)+'\t'+str(tnrate)+'\t'+str(fnrate)+'\t'+str(accuracy)+'\t'+str(taxamasked)+'\t'+str(percentmasked)+'\t'+str(percentgaps)+'\t'+case+'\t'+gene+'\tlinsi\tfubar\tzero\n')	
					
		for case in prefix:
			if case=='guidance' or case=='BMweights' or case=='PDweights':
				mask='70'
			else:
				mask='50'
				
			name = case+mask+'_'+str(n)+'.fasta'
				
			casealn=alndir+name
			casefu=fudir+name+'.fubar'
			parsed=AlignIO.read(casealn, 'fasta')
			
			testprobs = parseFubar(mapRef, casefu)	
				
			## Ensure that the mapping went ok.
			if len(truepos)!=len(testprobs):
				print case, len(map), len(truepos), len(testprobs)
				print "There is a serious flaw in your logic. Go cry in the corner.\n\n"
				assert 1==0

			(tp,tn,fp,fn,tprate,fprate,tnrate,fnrate,accuracy)=sweepRates(float(genes[gene]), truepos, testprobs)
			(taxamasked, percentmasked, percentgaps)=percentMaskedRes(parsed, alnlen, numseq)
			outhandle.write(str(n)+'\t'+str(tprate)+'\t'+str(fprate)+'\t'+str(tnrate)+'\t'+str(fnrate)+'\t'+str(accuracy)+'\t'+str(taxamasked)+'\t'+str(percentmasked)+'\t'+str(percentgaps)+'\t'+case+'\t'+gene+'\tlinsi\tfubar\t'+masks[mask]+'\n')	
outhandle.close()

			
			
			
		
				
			
			
			
			
			
			