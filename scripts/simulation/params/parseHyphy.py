#### SJS 8/9/13 ####

######### INCREDIBLY IMPORTANT NOTE: BEFORE RUNNING THIS SCRIPT ON A RAW HYPHY OUTPUT FILE, YOU'LL NEED TO DO A GREP FIND/REPLACE WITH THESE THREE FINDS: { } ^\s+ . REPLACE WITH NOTHING ###############

## Given a processed (SEE ABOVE COMMENT!!) output file from hyphy, run with the batch file REL_5cat.bf, determines a dN/dS estimate for each site.
# USAGE: python parseHyPhy.py <path_to_hyphy_output_file> <number_of_REL_categories> <outfile>

import re
import sys

#################################################################
def getSiteRates(post_matrix,  omegas):
	'''Create a list of rates at each position. Compute each as a weighted average = sum(postprob * w).'''

	allrates=[]
	for n in range(len(post_matrix[0])):
		rate = 0
		for w in range(0,len(omegas)):
			temp = float(postmat[w][n]) * float(omegas[w])	
			rate = rate + temp
		allrates.append(rate)
	return allrates
#################################################################

file = sys.argv[1]
ncat = int(sys.argv[2])
outfile = sys.argv[3]

assert(os.path.exists(file)), "Provided file doesn't exist. Check path?"
assert(ncat>=1), "You have no rate categories."

### Read in file
file = open(file, 'rU')
filestring = file.read()
file.seek(0)
lines=file.readlines()
numlines=len(lines)
file.close()

### Parse file for omega values
omegas=[]
for w in range(numcat+1):
	find=re.search('w_'+str(w)+'=(\d\.*\d*)', filestring)
	find_scinot = re.search('w_'+str(w)+'=\d\.*\d*e-\d+', filestring)
	if find != None:
		if find_scinot:
			omegas.append(float(0))  ## If small enough to resort to scientific notation, can round it to 0.
		else:
			omegas.append(float(find.group(1)))
			
## Get the posterior probability matrix of postprobs for site rates. 
#This is the THIRD matrix in the file (as in, the last one at the bottom)
postmat=[]
for n in range(numlines-numcat, numlines):
	entry=''
	probs=[]
	probs_raw=lines[n]
	for x in range(0, len(probs_raw)):
		if probs_raw[x]!=',' and probs_raw[x]!='\n':
			entry=entry+probs_raw[x]
		else:
			probs.append(float(entry))
			entry=''
	postmat.append(probs)

## Compute rates at each site by weighting omegas with post probs
rates = getSiteRates(postmat, omegas)

## Save to file
outhandle = open(outfile, 'w')
outhandle.write('site', 'omega')
sitecount=1
for r in rates:
	outhandle.write(str(sitecount)+'\t'+str(r)+'\n')
	sitecount+=1
outhandle.close()

	