# SJS 4/30/14.
# Function to count gaps in the reference alignments. Reported as percentages.

from Bio import AlignIO
import numpy as np

genes = ['or5', 'rho', 'prk', 'flat']
alndir = '/Users/sjspielman/Dropbox/aln/results/HA/alntree/aaguided'

for gene in genes:
	perc = np.zeros(100, dtype='float')
	dir = alndir + '_' + gene + '/'
	for i in range(100):
		aln = dir + "refaln" + str(i) + ".fasta"
		parsed = AlignIO.read(aln, 'fasta')
		for record in parsed:
			seq = str(record.seq)
			perc[i] = float(seq.count('-')) / float(len(seq))
	print np.mean(perc), np.std(perc)

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		