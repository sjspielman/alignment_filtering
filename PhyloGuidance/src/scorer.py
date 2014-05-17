from numpy import *
from Bio import AlignIO
import subprocess

class PrepScorer:
	def __init__(self):
		'''initialization function'''
		return	

	def countGaps(self, parsed, alnlen): #parsed = a parsed alignment.
		## Count the number of gaps in a column and returns a list. Does not give any information as to which positions are gaps. 		
		numGaps = []
		for i in range(0,alnlen): #number of sequences			
			findgaps = str(parsed[:,i])
			num = findgaps.count('-')
			numGaps.append(num)
		return numGaps
	
	def notGapSites(self, parsed, alnlen):
		## Determine which sites in a column are gaps. Returns a list of lists whereby each nested list contains the indices for which ones are not gaps or a 'G' for if that position is a gap
		notGaps=[]
		for i in range(0,alnlen): #number of sequences			
			findgaps = str(parsed[:,i])
			shortlist=[]
			for counter in range(len(findgaps)): #column!
				if findgaps[counter]!='-':
					shortlist.append(counter)
				else:
					shortlist.append('G')
			notGaps.append(shortlist)
		return notGaps
		
				
		
class ScoreProcessor(PrepScorer):
	def __init__(self):
		'''initialization function'''
		return

	
	
	###############################################################################################################
	def processScoresBMweights(self, n, numseq, alnlen, refMSA_file, ordered_weights, all_scores, allscores_file):
		
		#Parse reference alignment and count gaps per column. Save indices for which sequences ARE NOT GAPS, in notGaps
		parsed = AlignIO.read(refMSA_file, 'fasta')
		notGaps=self.notGapSites(parsed, alnlen) # Making a bunch of lists containing indices for if it's not a gap, 'G' if it is a gap. 
	
		# getting a normalization for each position. gap norms are 'G' (need placeholders)
		normScores=[] ##list of lists, for each column, of the normalization score at each position.
		for column in notGaps:
			shortlist=[]
			
			# Singletons first
			findgaps = column.count('G')
			if findgaps == len(column) - 1:
				[shortlist.append(0) for i in range(len(column))]
			
			else:
				for entry in column: # **Do not exclude self**
					if entry=='G':
						shortlist.append('G')
					else:
						thisweight=ordered_weights[entry] # get weight for position of interest
						norm=0
						for position in column:
							if position=='G': #skip the gaps
								continue
							else:
								norm=norm+(thisweight*ordered_weights[position])		
						shortlist.append(norm*n)
			normScores.append(shortlist)
		
		#Normalize scores. yes, this is gross.
		final_scores = zeros(numseq)
		col = 0
		for column in all_scores.T: #transposes it so loops over columns of all_scores instead of default rows
			colscores=zeros(numseq)
			taxon=0
			for rawscore in column:
				if normScores[col][taxon]=='G' or normScores[col][taxon] == 0:
					colscores[taxon]=0
				else:
					score=rawscore/(normScores[col][taxon])
					colscores[taxon]=score
				taxon+=1
			final_scores = vstack((final_scores, colscores)) #stack columnwise
			col+=1
		final_scores=final_scores.T
		final_scores = delete(final_scores, 0, axis=1)
		
		index_row=arange(1,(alnlen+1),1)
		final_scores_output = vstack((index_row, final_scores))
		savetxt(allscores_file, final_scores_output, delimiter='\t', fmt='%.3f')
		
		return final_scores	
	###############################################################################################################
	
		
		
	###############################################################################################################		
	def processScoresPDweights(self, n, numseq, alnlen, refMSA_file, dist_matrix, all_scores, allscores_file):

		#Parse reference alignment and count gaps per column. Save indices for which sequences ARE NOT GAPS, in notGaps
		parsed = AlignIO.read(refMSA_file, 'fasta')
		notGaps=self.notGapSites(parsed, alnlen) # Making a bunch of lists containing indices for if it's not a gap, 'G' if it is a gap. 
		
		# Contains the sum of pairwise distances for each row as a reference. 
		normSums=[]
		for s in range(numseq):
			normSums.append(sum(dist_matrix[s]))
		# getting a normalization for each position. gap norms are 'G' (need placeholders)
		# start with the sum of pairwise dists. then subtract out any comparisons to gaps.
		normScores=[] ## gives the normalization score for each column
		for column in notGaps:
			shortlist=[]
			
			# Deal with singletons.
			findgaps = column.count('G')
			if findgaps == len(column) - 1:
				[shortlist.append(0) for i in range(len(column))]
			
			else:
				for i in column: # i is also the index.
					if i == 'G':
						shortlist.append('G')
					else:
						newNorm = normSums[i]
						count=0
						for x in column:
							if x =='G':
								newNorm = newNorm - dist_matrix[count][i]
							count+=1
						shortlist.append(newNorm*n)	
			normScores.append(shortlist)
		#Normalize scores.
		final_scores = zeros(numseq)
		col=0
		for column in all_scores.T: #transposes it so loops over columns of all_scores instead of default rows
			colscores=zeros(numseq)
			taxon=0
			for rawscore in column:
				if normScores[col][taxon]=='G' or normScores[col][taxon] == 0: # give scores of 0 to gaps or to singletons (which are assigned a normalization of 0)
					colscores[taxon]=0
				else:
					score=rawscore/(normScores[col][taxon])
					colscores[taxon]=score
				taxon+=1
			final_scores = vstack((final_scores, colscores)) #stack columnwise
			col+=1
		final_scores=final_scores.T
		final_scores = delete(final_scores, 0, axis=1)
		
		index_row=arange(1,(alnlen+1),1)
		final_scores_output = vstack((index_row, final_scores))
		savetxt(allscores_file, final_scores_output, delimiter='\t', fmt='%.3f')
	
		return final_scores
	###############################################################################################################
		
	
	##############################################################################################################
	def processScoresGuidance(self, n, numseq, alnlen, refMSA_file, all_scores, allscores_file):
	
		parsed = AlignIO.read(refMSA_file, 'fasta')
		numGaps=self.countGaps(parsed, alnlen)
	
		final_scores = zeros(alnlen)

		#Normalize scores
		normedscores=zeros(alnlen)
		t=0 #taxon counter
		s=0 #entry in row counter
		for taxon in all_scores:
			for score in taxon:
				norm = (numseq - numGaps[s] -1)*n   ## the "-1" is because we can't be comparing it to itself.
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
		
		index_row=arange(1,(alnlen+1),1)
		final_scores_output = vstack((index_row, final_scores))
		savetxt(allscores_file, final_scores_output, delimiter='\t', fmt='%.3f')
		
		return (final_scores)
	###############################################################################################################

	### ALGORITHMS WHICH PENALIZE GAPS
	
	def processScoresPenalize(self, n, numseq, alnlen, all_scores, allscores_file):
		normScores = n*(numseq-1)
		final_scores = all_scores/normScores
		
		index_row=arange(1,(alnlen+1),1)
		final_scores_output = vstack((index_row, final_scores))
		savetxt(allscores_file, final_scores_output, delimiter='\t', fmt='%.3f')
		
		return final_scores
	
	
	def processScoresPenalizeBMweights(self, n, alnlen, ordered_weights, all_scores, allscores_file):
		final_scores = zeros(alnlen)
		normedscores = zeros(alnlen)
		
		count=0
		for taxon in all_scores:
			normedscores = taxon/(ordered_weights[count] * n)
			final_scores = vstack((final_scores, normedscores))
			normedscores=zeros(alnlen)
			count+=1
		final_scores = delete(final_scores, (0), axis=0)
		
		index_row=arange(1,(alnlen+1),1)
		final_scores_output = vstack((index_row, final_scores))
		savetxt(allscores_file, final_scores_output, delimiter='\t', fmt='%.3f')
		
		return final_scores
	
	def processScoresPenalizePDweights(self, n, alnlen, dist_matrix, all_scores, allscores_file):
		final_scores = zeros(alnlen)
		normedscores = zeros(alnlen)
				
		count=0
		for taxon in all_scores:
			norm = float(sum(dist_matrix[count]))
			normedscores = taxon/(norm * n)
			final_scores = vstack((final_scores, normedscores))
			normedscores=zeros(alnlen)
			count+=1
		final_scores = delete(final_scores, (0), axis=0)
		
		index_row=arange(1,(alnlen+1),1)
		final_scores_output = vstack((index_row, final_scores))
		savetxt(allscores_file, final_scores_output, delimiter='\t', fmt='%.3f')
		
		return final_scores





class Scorer(PrepScorer,ScoreProcessor):
	def __init__(self):
		'''initialization function'''
		return

	def scoreMSA_Guidance(self, refMSA_file, n, numseq, alnlen, allscores_file, allscores_file_penalized, numSaveTrees):
		n_trees= len(numSaveTrees)
		alg='g'
		all_scores=zeros((numseq, alnlen))		
		for i in range(n_trees):
			testMSA_file='bootaln'+str(i)+'.fasta'
			outfile = "scores"+str(i)+"_"+str(alg)+".txt"
			scoreCommand='../src/score/guidance_score ' + refMSA_file + " " + testMSA_file + " " + outfile
			subprocess.call(scoreCommand, shell=True)
			scores = loadtxt(outfile)
			
			# Multiply the scores by how many times numSaveTrees calls for
			scores = scores*numSaveTrees[i]	
			all_scores = all_scores + scores
		
		final_scores=self.processScoresGuidance(n, numseq, alnlen, refMSA_file, all_scores, allscores_file)
		final_scores_penalized=self.processScoresPenalize(n, numseq, alnlen, all_scores, allscores_file_penalized)
		return (final_scores, final_scores_penalized)
	
	
	def scoreMSA_BMweights(self, refMSA_file, n, numseq, alnlen, ordered_weights, weightfile, allscores_file, allscores_file_penalized, numSaveTrees):
		n_trees = len(numSaveTrees)
		alg='bm'
		all_scores=zeros((numseq, alnlen))
		for i in range(n_trees):
			testMSA_file='bootaln'+str(i)+'.fasta'
			outfile = "scores"+str(i)+"_"+str(alg)+".txt"
			scoreCommand='../src/score/bmweights_score ' + refMSA_file + " " + testMSA_file + " " + weightfile + " " + outfile
			subprocess.call(scoreCommand, shell=True)
			scores = loadtxt(outfile)
			
			# Multiply the scores by how many times numSaveTrees calls for
			scores = scores*numSaveTrees[i]	
			all_scores = all_scores + scores
			
		final_scores=self.processScoresBMweights(n, numseq, alnlen, refMSA_file, ordered_weights, all_scores, allscores_file)
		final_scores_penalized = self.processScoresPenalizeBMweights(n, alnlen, ordered_weights, all_scores, allscores_file_penalized)
		return (final_scores, final_scores_penalized)
	
	
	def scoreMSA_PDweights(self, refMSA_file, n, numseq, alnlen, dist_matrix, dist_matrix_file, allscores_file, allscores_file_penalized, numSaveTrees):
		n_trees = len(numSaveTrees)
		alg='pd'
		all_scores=zeros((numseq, alnlen))
		for i in range(n_trees):
			testMSA_file='bootaln'+str(i)+'.fasta'
			outfile = "scores"+str(i)+"_"+str(alg)+".txt"
			scoreCommand='../src/score/pdweights_score ' + refMSA_file + " " + testMSA_file + " " + dist_matrix_file + " " + outfile
			subprocess.call(scoreCommand, shell=True)
			scores = loadtxt(outfile)
			
			# Multiply the scores by how many times numSaveTrees calls for
			scores = scores*numSaveTrees[i]					
			all_scores = all_scores + scores
			
		final_scores = self.processScoresPDweights(n, numseq, alnlen, refMSA_file, dist_matrix, all_scores, allscores_file)
		final_scores_penalized = self.processScoresPenalizePDweights(n, alnlen, dist_matrix, all_scores, allscores_file_penalized)
		return (final_scores, final_scores_penalized)
	