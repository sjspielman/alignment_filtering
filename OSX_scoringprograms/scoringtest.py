#!/usr/bin/python

from numpy import *
from Bio import SeqIO, AlignIO

def nogapPositions( seq ):
	'''This function returns an iterator that iterates over all positions in a sequence that don't correspond to a gap (-)'''
	pos = -1
	for c in seq:
		pos = pos + 1
		if c=='-': pass
		else: yield pos

def calcScores( base_MSA, test_MSA ):
	'''This function calculates a matrix of scores for the base MSA, relative to the test MSA. Caution: no sanity checking is performed. If the test MSA does not contain the exact same sequences as the base MSA, the function will fail.'''


	bc = len( base_MSA[1] ) # number of columns in base alignment
	tc = len( test_MSA[1] ) # number of columns in test alignment
	n = len( base_MSA ) # number of sequences

	# set up the map of positions from test to base array
	# by default, gaps map to -1
	test_to_base_map = -1*ones( (n, tc), dtype=int )

	ts_iter = test_MSA.__iter__() # iterator over test MSA sequences

	j = 0 # counts current sequence
	for bseq in base_MSA: # loop over all sequences in base MSA
		tseq = ts_iter.next()

		# loop over all no-gap positions in sequence from base MSA
		# and map to corresponding positions in test MSA
		ti_iter = nogapPositions( tseq )
		for bi in nogapPositions( bseq ):
			ti = ti_iter.next()
			test_to_base_map[ j, ti ] = bi
		j += 1

	print "\nMap from test MSA to base MSA:"
	print test_to_base_map # completed map

	# now that we have the map, we can calculate the scores in the original MSA
	scores = zeros( (n, bc) ) # start with matrix of zeros
	for ti in range( tc ): # go over all columns in map
		#print "Column", ti
		for j in range( n ): # go over all rows
			#print "Row", j
			x0 = test_to_base_map[ j, ti ]
			# sum up all the pairs that point to the same location in base MSA
			# we have to subtract 1 because each site is also compared to itself
			if x0 >= 0:
				scores[ j, x0 ] += sum([1 for x in test_to_base_map[ :, ti ] if x == x0 ]) - 1	
	return scores

def calcScores2( base_MSA, test_MSA ):
	'''This function calculates a matrix of scores for the base MSA, relative to the test MSA. Caution: no sanity checking is performed. If the test MSA does not contain the exact same sequences as the base MSA, the function will fail.'''


	bc = len( base_MSA[1] ) # number of columns in base alignment
	tc = len( test_MSA[1] ) # number of columns in test alignment
	n = len( base_MSA ) # number of sequences

	# set up the map of positions from test to base array
	# by default, gaps map to -1
	test_to_base_map = -1*ones( (n, tc), dtype=int )

	ts_iter = test_MSA.__iter__() # iterator over test MSA sequences

	j = 0 # counts current sequence
	for bseq in base_MSA: # loop over all sequences in base MSA
		tseq = ts_iter.next()

		# loop over all no-gap positions in sequence from base MSA
		# and map to corresponding positions in test MSA
		ti_iter = nogapPositions( tseq )
		for bi in nogapPositions( bseq ):
			ti = ti_iter.next()
			test_to_base_map[ j, ti ] = bi
		j += 1

	print "\nMap from test MSA to base MSA:"
	print test_to_base_map # completed map

	# now that we have the map, we can calculate the scores in the original MSA
	scores = zeros( (n, bc) ) # start with matrix of zeros
	for ti in range( tc ): # go over all columns in map
		#print "Column", ti
		for j in range( n ): # go over all rows
			#print "Row", j
			x0 = test_to_base_map[ j, ti ]
			# sum up all the pairs that point to the same location in base MSA
			# we have to subtract 1 because each site is also compared to itself
			if x0 >= 0:
				scores[ j, x0 ] += sum([test_to_base_map[ :, ti ] == x0 ]) - 1	
	return scores

def calcScores3( base_MSA, test_MSA ):
	'''This function calculates a matrix of scores for the base MSA, relative to the test MSA. Caution: no sanity checking is performed. If the test MSA does not contain the exact same sequences as the base MSA, the function will fail.'''


	bc = len( base_MSA[1] ) # number of columns in base alignment
	tc = len( test_MSA[1] ) # number of columns in test alignment
	n = len( base_MSA ) # number of sequences

	# set up the map of positions from test to base array
	# by default, gaps map to bc (one past the last column in the base alignment)
	test_to_base_map = bc*ones( (n, tc), dtype=int )

	ts_iter = test_MSA.__iter__() # iterator over test MSA sequences

	j = 0 # counts current sequence
	for bseq in base_MSA: # loop over all sequences in base MSA
		tseq = ts_iter.next()

		# loop over all no-gap positions in sequence from base MSA
		# and map to corresponding positions in test MSA
		ti_iter = nogapPositions( tseq )
		for bi in nogapPositions( bseq ):
			ti = ti_iter.next()
			test_to_base_map[ j, ti ] = bi
		j += 1

	print "\nMap from test MSA to base MSA:"
	print test_to_base_map # completed map

	# now that we have the map, we can calculate the scores in the original MSA
	scores = zeros( (n, bc+1) ) # start with matrix of zeros
	# we add an extra column to store scores for gaps, that colum will be deleted at the end
	for j in range( n ):  # go over all rows
		for ti in range( tc ): # go over all columns in map
			x0 = test_to_base_map[ j, ti ]
			#print x0
			# sum up all the pairs that point to the same location in base MSA
			# we have to subtract 1 because each site is also compared to itself
			scores[ j, x0 ] += sum([test_to_base_map[ :, ti ] == x0 ]) - 1
	return scores[ :, 0:bc]

def calcScores4( base_MSA, test_MSA ):
	'''This function calculates a matrix of scores for the base MSA, relative to the test MSA. Caution: no sanity checking is performed. If the test MSA does not contain the exact same sequences as the base MSA, the function will fail.'''


	bc = len( base_MSA[1] ) # number of columns in base alignment
	tc = len( test_MSA[1] ) # number of columns in test alignment
	n = len( base_MSA ) # number of sequences

	# set up the map of positions from test to base array
	# by default, gaps map to bc (one past the last column in the base alignment)
	test_to_base_map = bc*ones( (n, tc), dtype=int )

	ts_iter = test_MSA.__iter__() # iterator over test MSA sequences

	j = 0 # counts current sequence
	for bseq in base_MSA: # loop over all sequences in base MSA
		tseq = ts_iter.next()

		# loop over all no-gap positions in sequence from base MSA
		# and map to corresponding positions in test MSA
		ti_iter = nogapPositions( tseq )
		for bi in nogapPositions( bseq ):
			ti = ti_iter.next()
			test_to_base_map[ j, ti ] = bi
		j += 1

	print "\nMap from test MSA to base MSA:"
	print test_to_base_map # completed map

	# now that we have the map, we can calculate the scores in the original MSA
	scores = zeros( (n, bc+1) ) # start with matrix of zeros
	# we add an extra column to store scores for gaps, that column will be deleted at the end
	for j in range( n ):  # go over all rows
		x0 = test_to_base_map[ j, : ] # vector of indices
		scores[ j, x0] = [sum([test_to_base_map[ :, i ] == x0[i] ]) - 1 for i in range(tc)]
	return scores[ :, 0:bc]
	
def quicktest():
	base_MSA = ['ACCTG','-CCTG','ACCAG','-CCA-']
	test_MSA = ['A-CCTG','--CCTG','A-CCAG','-C-CA-']

	print "Base MSA:"
	for seq in base_MSA: print "   ", seq

	print "\nTest MSA:"
	for seq in test_MSA: print "   ", seq

	scores = calcScores( base_MSA, test_MSA )
	print "\nMatrix of scores:"
	print scores
	
	scores = calcScores4( base_MSA, test_MSA )
	print "\nMatrix of scores:"
	for row in scores:
		for i in row:
			print i,
		print


def longtest():
	handle = open("full.aln", "rU")
	base_MSA = []
	for record in SeqIO.parse(handle, "fasta") :
		base_MSA.append( record.seq )
	handle.close()
	
	test_MSA = base_MSA
	scores = calcScores4( base_MSA, test_MSA )
	print "\nMatrix of scores:"
	for row in scores:
		for i in row:
			print i,
		print

msa=AlignIO.read('BP_MSA/MSA.MAFFT.aln', 'fasta')
baseMSA=[]
for entry in msa:
	baseMSA.append(entry.seq)
msa=AlignIO.read('BP_MSA/MSA.MAFFT.tree_1.MAFFT.aln', 'fasta')
testMSA=[]
for entry in msa:
	testMSA.append(entry.seq)
print baseMSA

#quicktest()
allscores_file='calcScores4_run1.txt'
scores=calcScores( baseMSA, testMSA )
#savetxt(allscores_file, scores, delimiter='\t', fmt='%.6f')
