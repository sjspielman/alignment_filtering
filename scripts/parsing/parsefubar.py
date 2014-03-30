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
	
	# Now find the true evolutionary rates using mapTrue
	truepos =  parseTrueRates(trfile, truealn, mapTrue)

	return (mapRef, truepos)
	
	
def getConsensus(allMaps):
	consMap=[] # will have final consensus map
	
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
################################################################################################
################################################################################################
################################################################################################
################################################################################################


prefix=['refaln', 'Guidance', 'BMweights', 'PDweights', 'GuidanceP', 'BMweightsP', 'PDweightsP']
genes=['rho']
mask={'50':'fifty'}

datadir='/Users/sjspielman/Dropbox/aln/results/'
if type=='neutral':
	datadir += 'neutral/'
elif type == 'HA':
	datadir += 'HA/'


outfile='fubar90.txt'
outhandle=open(outfile, 'w')
outhandle.write('count\ttprate\tfprate\t\tfnrate\taccuracy\tcase\tgene\tmask\tmethod\tpenal\n')


for gene in genes:


	print gene+'\n'
	
	# Directories: fubar output, alignments (all made with linsi)
	fudir = datadir+'fubar/fubar_'+gene+'/'
	alndir = datadir+'alntree/nucguided_'+gene+'/'
	
	# Directories: true simulated alignments and evolutionary rate categories
	truerates_dir=datadir+'Simulation/truerates/'+gene+'/'
	truealn_dir=datadir+'Simulation/sequences/'+gene+'/'
	
			
	for n in range(100):
		print str(n)
		
		# True alignment and evolutionary rate files
		trfile=truerates_dir+'truerates'+str(n)+'.txt'
		truealn=truealn_dir+'truealn_codon'+str(n)+'.fasta'
		trueparsed=AlignIO.read(truealn, 'fasta')
		
		# Reference alignment
		refaln=alndir+'refaln'+str(n)+'.fasta'
		refparsed=AlignIO.read(refaln, 'fasta')
		alnlen=len(refparsed[0])
		numseq=len(refparsed)
			
		# Build map to true alignment and obtain simulated positive selection state (binary - 0=notpos, 1=pos)
		(mapRef, truepos) = buildMap(trueparsed, refparsed, numseq, alnlen)
		
		
		########### Accuracy assessment #########
		for mask in masks:
			for case in prefix:
				
				## Get file names and whether or not gap-penalized algorithm
				if case=='refaln':
					penal='zero'
					fubar=fudir+'refaln'+str(n)+'.fasta.fubar'
					aln=refaln
					parsed=refparsed
				
				elif case=='Guidance' or case=='BMweights' or case=='PDweights':
					penal='no'
					name = case+mask+'_'+str(n)+'.fasta'
					aln=alndir+name		
					fubar=fudir+name+'.fubar'
					parsed=AlignIO.read(aln, 'fasta')
				
				else:
					penal='yes'
					name = case+mask+'_'+str(n)+'.fasta'
					aln=alndir+name
					fubar=fudir+name+'.fubar'
					parsed=AlignIO.read(aln, 'fasta')	
					
				testprobs = parseFubar(mapRef, fubar)	
						
				## Ensure that the mapping went ok.
				if len(truepos)!=len(testprobs):
					print case, len(map), len(truepos), len(testprobs)
					print "Mapping failed."
					assert 1==0
	
	
				## Fubar assessment	at single posterior probability cutoff			
				(tp,tn,fp,fn,tprate,fprate,tnrate,fnrate,accuracy)=sweepRates(0.895, truepos, testprobs)
				outhandle.write(str(n)+'\t'+str(tprate)+'\t'+str(fprate)+'\t'+str(fnrate)+'\t'+str(accuracy)+'\t'+case+'\t'+gene+'\t'+masks[mask]+'\tfubar\t'+penal+'\n')	

outhandle.close()

			
			
			
		
				
			
			
			
			
			
			