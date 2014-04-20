import re, sys, csv
from numpy import *
from Bio import AlignIO, SeqIO

### SJS ###
# Contains all the functions used by the fubar/paml parsing/accuracy inference scripts.


###########################################################################################################
def singleTaxonMap(trueparsed, parsed, numseq, alnlen):
	''' Map for relevant sites for both the true (given by Indelible) and reference (inferred, no filtering) based on a single taxon. '''
	# trueparsed: the true alignment
	# parsed: the reference alignment
	# numseq: number of sequences in alignment.
	# alnlen: alignment length of the REFERENCE alignment.
	
	reftaxon = '11' # Arbitrary, based entirely on my birthday.
	
	# Get the true and reference alignment sequences for this taxon
	trueseq = ''
	refseq = ''
	for entry in trueparsed:
		if str(entry.id)==reftaxon:
			trueseq=str(entry.seq)
			break
	for entry in parsed:
		if str(entry.id)==reftaxon:
			refseq=str(entry.seq)
			break
	#Build the map. Record the index for each non-gap site in trueseq. Then, go through the refseq and for the nongaps, pop off that index. For the gaps, add a '-1'.
	truelist=[]
	map2True = []
	for a in range(len(trueseq)):
		if trueseq[a] != '-': 
			truelist.append(a)	### Contains positions which are NOT gaps in the true alignment	
	#print truelist
	for a in range(len(refseq)):
		if refseq[a] == '-':
			map2True.append(-1)    
		else:
			map2True.append(truelist.pop(0))
	#print map2True
	
	# Convert.
	wantRef = []
	wantTrue = []
	for i in range(len(map2True)):
		# if refseq is a gap, we don't want 
		if map2True[i] != -1:
			wantRef.append(i)
			wantTrue.append(map2True[i])
	assert ( len(wantRef) > 0 and len(wantRef) == len(wantTrue) ), "Bad map building."
	return wantRef, wantTrue
###########################################################################################################

###########################################################################################################
def consensusMap(trueparsed, parsed, numseq, alnlen):
	''' Build maps to relevant sites for both the true (given by Indelible) and reference (inferred, no filtering) based on >=50% of sites being in that column. '''
	# trueparsed: the true alignment
	# parsed: the reference alignment
	# numseq: number of sequences in alignment.
	# alnlen: alignment length of the REFERENCE alignment.
	
	allMaps = zeros(alnlen, dtype=int)
	
	for entry in range(numseq):
	
		taxon = str(entry) ### Taxon names are integers 0-numseq
		# Will contain indices for positions in the TRUE alignment sequence which are NOT gaps
		truelist=[]                        
		
		## The index in map2True corresponds to the refaln position. The value in map2True corresponds to the truealn position. If the reference position is a gap, the value will instead be -1.
		map2True=zeros(alnlen, dtype=int)
			
		# Get the true and reference alignment sequences for this taxon
		trueseq = ''
		refseq = ''
		for entry in trueparsed:
			if str(entry.id)==taxon:
				trueseq=str(entry.seq)
				break
		for entry in parsed:
			if str(entry.id)==taxon:
				refseq=str(entry.seq)
				break
		
		#Build the map. Record the index for each non-gap site in trueseq. Then, go through the refseq and for the nongaps, pop off that index. For the gaps, add a '-1'.
		for a in range(len(trueseq)):
			if trueseq[a] != '-': 
				truelist.append(a)		
		for a in range(len(refseq)):
			if refseq[a] == '-':
				map2True[a] = '-1'    
			else:
				map2True[a]=truelist.pop(0)
		
		allMaps = vstack((allMaps, map2True))
	allMaps = delete(allMaps, 0, axis=0)	

	# Build a 50% consensus map.
	counter = 0
	wantRef = []
	wantTrue = []
	for column in allMaps.T: # transpose so can loop over columns.
		
		# Check that the most frequent column entry occurs at least 50% of the time.			
		col = column.tolist()
		most = int(max(set(col), key=col.count))
		countmost = col.count(most)
		perc_max = float(countmost)/float(len(col))
		if perc_max > 0.5: # Needs to be >, not >= to ensure columns aren't repeated
			wantRef.append(counter) # "counter" refers to the refaln position
			wantTrue.append(most)   # "most" refers to the truealn position 
		counter+=1
	#print float(len(wantTrue))/float(alnlen)
	return (wantRef, wantTrue)
###########################################################################################################

###########################################################################################################
def assessTrueFUBAR(trfile, fufile, posStart):
	''' This way no need to get tripped up by maps '''
	# trfile: file with true simulated rates, given by Indelible
	# posStart: omega category in which positive selection begins. Depends on simulation.
	poslist=[]
	
	infile=open(trfile, 'r')
	truelines=infile.readlines()
	infile.close()
	truelines=truelines[10:] ## only keep these lines since before that it's all header crap.	

	for counter in range(len(truelines)):
		find=re.search('^\d+\t(\d+)\t', truelines[counter])
		assert(find), "Could not parse truerates file."
		if find:
			rate=int(find.group(1))
			if rate >=posStart: 
				poslist.append(1)
			else:
				poslist.append(0)
	
	## Read in fubar file
	testprobs=[] ## contains the prob(alpha>beta) values for the truefubar results
	fubar=csv.reader(open(fufile,'r'))
	allfubar=[]
	
	for row in fubar:
		if row[0]=='Codon':
			continue
		else:
			allfubar.append(float(row[4]))
	
	# corresponds to truepos, testprobs
	return(poslist, allfubar)
###########################################################################################################
	
###########################################################################################################
def assessTruePAML(trfile, fufile, posStart):
	''' This way no need to get tripped up by maps '''
	# trfile: file with true simulated rates, given by Indelible
	# posStart: omega category in which positive selection begins. Depends on simulation.
	poslist=[]
	
	infile=open(trfile, 'r')
	truelines=infile.readlines()
	infile.close()
	truelines=truelines[10:] ## only keep these lines since before that it's all header crap.	

	for counter in range(len(truelines)):
		find=re.search('^\d+\t(\d+)\t', truelines[counter])
		assert(find), "Could not parse truerates file."
		if find:
			rate=int(find.group(1))
			if rate >=posStart: 
				poslist.append(1)
			else:
				poslist.append(0)
	
	# read paml file
	allpaml=[]
	paml = open(paml_file, 'r')
	all_lines = paml.readlines()
	paml.close()
	
	prior, omega = getMLE(all_lines)
	
	
	#First par down the PAML file to keep only the relevant lines
	counter=0
	for line in all_lines:
		find = re.search('^Bayes Empirical Bayes', line)
		if find:
			lines = all_lines[counter+3:counter+3+alnlen]
			break
		else:
			counter+=1
			continue
	
	# Retrieve pr(w>1)
	for line in lines:
		find = re.search('^\s*\d+\s.+ (\d\.\d+) \(\s*\d+\)', line)
		if find:
			# just to kill scientific notation. if it's this small who cares what it's prob is. doesn't even matter for sweeping.
			if float(find.group(1)) <= 0.001:
				allpaml.append(0.001)
			else:
				allpaml.append(float(find.group(1)))
	
	# corresponds to truepos, testprobs
	return(poslist, allpaml, prior, omega)
###########################################################################################################

###########################################################################################################
def parsePAML(map, paml_file, alnlen):
	allpaml=[]
	testprobs=[]
	
	paml = open(paml_file, 'r')
	all_lines = paml.readlines()
	paml.close()
	
	prior, omega = getMLE(all_lines)
	
	#First par down the PAML file to keep only the relevant lines
	counter=0
	for line in all_lines:
		find = re.search('^Bayes Empirical Bayes', line)
		if find:
			lines = all_lines[counter+3:counter+3+alnlen]
			break
		else:
			counter+=1
			continue
	
	# Retrieve pr(w>1)
	for line in lines:
		find = re.search('^\s*\d+\s.+ (\d\.\d+) \(\s*\d+\)', line)
		if find:
		
			# just to kill scientific notation. if it's this small who cares what it's prob is. doesn't even matter for sweeping.
			if float(find.group(1)) <= 0.001:
				allpaml.append(0.001)
			else:
				allpaml.append(float(find.group(1)))
	for entry in map:
		testprobs.append(allpaml[entry])	
	return (testprobs, prior, omega)
###########################################################################################################

###########################################################################################################
def getMLE(paml):
	''' Get the prior (line 6) and omega (line 7) for the w>1 category ''' 
	#pline = paml[5]
	#wline = paml[6]
	
	findp = re.search("\s+(\d\.\d+)$", paml[5])
	assert (findp), "Could not find the prior probability from PAML."
	prior = findp.group(1)
	
	findw = re.search("\s+(\d+\.\d+)$", paml[6])
	assert (findw), "Could not find the prior probability from PAML."
	omega = findw.group(1)

	return (prior, omega)
###########################################################################################################

###########################################################################################################	
def gridWeights(fubar):
	''' Something with grid weights. '''
	hi = 1
	return 1
###########################################################################################################

###########################################################################################################
def parseFUBAR(map, fufile):
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
###########################################################################################################
		
###########################################################################################################
def parseTrueRates(trfile, wantTrue, posStart):
	''' Retrieve binary list (0=not pos, 1=pos) for true simulated rates at sites of interest ''' 
	# trfile: file with true simulated rates, given by Indelible
	# wantTrue: relevant sites in true alignment
	# posStart: omega category in which positive selection begins. Depends on simulation.

	# For positions we care about, 0 if omega<=1 , 1 if omega>1
	poslist=[]
	
	infile=open(trfile, 'r')
	truelines=infile.readlines()
	infile.close()
	truelines=truelines[10:] ## only keep these lines since before that it's all header crap.	

	for counter in wantTrue:
		#print truelines[counter]
		find=re.search('^\d+\t(\d+)\t', truelines[counter])
		assert(find), "Could not parse truerates file."
		if find:
			rate=int(find.group(1))
			if rate >=posStart: 
				poslist.append(rate)
			else:
				poslist.append(0)
	return poslist
###########################################################################################################


###########################################################################################################
def getAccuracy(x, truepos, testprobs):
	'''Determine accuracy statistics (tpr, fpr, etc).'''
	# x: posterior probability cutoff
	# truepos: whether that position was positively selected or not in the true alignment
	# testprobs: inferred probability that the site is positively selected
	
	tp, tn, fn, fp = 0, 0, 0, 0
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
	(tprate,fprate,tnrate,fnrate,accuracy) = calcStats(tp, fp, tn, fn)

	return (tp,tn,fp,fn,tprate,fprate,tnrate,fnrate,accuracy)
###########################################################################################################

###########################################################################################################
def calcStats(tp, fp, tn, fn):
	''' See function name for details. '''
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
###########################################################################################################

###########################################################################################################
def assessMasking(alnfile):
	''' Determine extent to which alignment was filtering. Returns 3 values: num, ave, perc: '''
	## num     = total number of sites masked in the entire alignment
	## ave     = average number of sites masked per taxon
	## perc    = percent of residues (excludes gaps) which were masked in entire alignment
		
	parsed = AlignIO.read(alnfile, 'fasta')
	
	seq_string = ''
	seq_length = 0
	mask_per_taxon = zeros(len(parsed))
	perc_per_taxon = zeros(len(parsed))
	
	count = 0
	for record in parsed:
		seq = re.sub('-','',str(record.seq)) # remove gaps
		seq_string += seq
		seq_length += len(seq)
		mask_per_taxon[count] = seq.count('?') # save number of masks
		perc_per_taxon[count] = float(seq.count('?')) / float(len(seq)) # save length of this sequence	
		count += 1
	
	# Calc masking statistics
	num = seq_string.count('?')
	ave = mean(mask_per_taxon)
	perc = float(num) / float(seq_length)
	
	return (num, ave, perc)
###########################################################################################################


























###########################################################################################################
