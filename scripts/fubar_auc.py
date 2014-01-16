
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
		#Negative cases
		elif float(testprobs[i])<x:
			if truepos[i]==0:
				tn+=1
			elif truepos[i]>0:
				fn+=1
	(tprate,fprate) = calcStats(tp, fp, tn, fn)
	return (tprate,fprate)

def calcStats(tp, fp, tn, fn):
	'''Pretty obvious.'''
	if fn==0 and tp==0:
		#fnrate=0
		tprate=0
	else:
		#fnrate=(float(fn))/(float(fn+tp))
		tprate=(float(tp))/(float(fn+tp))
	if tn==0 and fp==0:
		#tnrate=0
		fprate=0
	else:
		#tnrate=(float(tn))/(float(tn+fp))
		fprate=(float(fp))/(float(tn+fp))
	#accuracy = (float(tp+tn)/float(tp+tn+fp+fn))
	return (tprate,fprate)

def calcAUC(truepos, testprobs, cutoffs):
	auc=0
	fprate_old=0			
	for x in cutoffs:
		(tprate, fprate)=sweepRates(x, truepos, testprobs)
	
		if fprate==1:
			auc += tprate_old * (limit-fprate_old)
			break
		
		if (fprate>limit and fprate_old!=0.01):
			auc += tprate_old * (limit-fprate_old)
			break
		
		auc += tprate * (fprate - fprate_old)
		
		if fprate==0.01:
			break
		
		fprate_old=fprate
		tprate_old=tprate
	return auc
#########################################################################################	
#########################################################################################	

prefix=['guidance', 'BMweights', 'PDweights', 'guidance_p', 'BMweights_p', 'PDweights_p']
genes=['or5','rho', 'prk', 'flat']
mask='50'
	

cutoffs=arange(1,-0.001,-0.001)
limit=0.01


outfile='analysis/fubar_auc.txt'
outhandle=open(outfile, 'w')
outhandle.write('count\tauc\tcase\tgene\tmethod\n')## count  auc  case  gene  method

for gene in genes:
	print gene

	fudir = 'fubar/fubar_'+gene+'_seqs_real/'
	alndir = 'alntree/nucguided_linsi_'+gene+'_seqs_real/'
	
	truerates_dir='../Simulation/truerates/'+gene+'/seqs_real/'
	truealn_dir='../Simulation/sequences/'+gene+'/seqs_real/'
	
	for n in range(100):
		print '\n'+str(n)
			
		trfile=truerates_dir+'truerates'+str(n)+'.txt'
		truealn=truealn_dir+'truealn_codon'+str(n)+'.fasta'
			
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
	
		## First process reference alignment
		testprobs = parseFubar(mapRef, reffu)
		
		## Ensure that the mapping went ok.
		if len(truepos)!=len(testprobs):
			print case, len(map), len(truepos), len(testprobs)
			print "There is a serious flaw in your logic. Go cry in the corner.\n\n"
			assert 1==0
		
		auc=calcAUC(truepos, testprobs, cutoffs)
		outhandle.write(str(n)+'\t'+str(auc)+'\trefaln\t'+gene+'\tfubar\tzero\n')	
		
		
		for case in prefix:
			print case
			if case=='guidance' or case=='BMweights' or case=='PDweights':
				penal='no'
			else:
				penal='yes'
				
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

			auc=calcAUC(truepos, testprobs, cutoffs)
			outhandle.write(str(n)+'\t'+str(auc)+'\t'+case+'\t'+gene+'\tfubar\t'+penal+'\n')	
outhandle.close()
			
			
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	