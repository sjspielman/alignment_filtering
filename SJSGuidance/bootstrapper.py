import subprocess, sys, shutil, os
from Bio import SeqIO
from numpy import *

class Bootstrapper(object): 
	def __init__(self, aligner, tree_builder, scorer):
		'''initialization function'''
		self.aligner = aligner
		self.tree_builder = tree_builder
		self.scorer = scorer
	
	def bootstrap( self, prealn, refMSA_file, n, numprocesses):
		
		# Obtain reference alignment, as well as number of sequences and length of alignment
		refaln_seq=[]
		infile=open(refMSA_file, 'r')
		parsed = list(SeqIO.parse(infile, 'fasta'))
		infile.close()
		for record in parsed:
			refaln_seq.append(str(record.seq))
		numseq = len(refaln_seq)
		alnlen = len(refaln_seq[0])	
		
		
		final_treefile = 'BStrees.tre'
		
		# Create the bootstrapped trees
		print "Constructing bootstrap trees"
		self.tree_builder.buildBootTrees(n, refaln_seq, numseq, alnlen, final_treefile)
			
		# Process trees for given alignment software
		print "Formatting trees"
		self.aligner.processTrees(n, final_treefile) 	
	
		# Create the boostrapped alignments, using the bootstrapped trees		
		print "Building bootstrap alignments"
		self.aligner.multiMakeAlignmentsGT(prealn, n, numprocesses)
		
		return(numseq, alnlen)
		


class guidanceBootstrapper(Bootstrapper):
	''' Score with original guidance only.'''
	
	def __init__(self, aligner, tree_builder, scorer):
		'''initialization function'''
		super(guidanceBootstrapper, self).__init__(aligner, tree_builder, scorer)

	def runBootstrap(self, BootDir, prealn, refMSA_file, pflag, n, numprocesses, allscores_file):
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn, refMSA_file, n, numprocesses)
		
		print "Scoring, original Guidance"
		gscores= self.scorer.scoreMSA_Guidance(refMSA_file, n, numseq, alnlen, pflag, allscores_file)

		return(numseq, alnlen, gscores)


class bmBootstrapper(Bootstrapper):
	''' Score with BranchManager weights only.'''
	
	def __init__(self, aligner, tree_builder, weight_tree_builder, scorer):
		'''initialization function'''
		self.weight_tree_builder = weight_tree_builder
		super(bmBootstrapper, self).__init__(aligner, tree_builder, scorer)
		
	def runBootstrap(self, BootDir, prealn, refMSA_file, pflag, n, numprocesses, allscores_file): 
		
		# Create the scoring tree
		weightfile='weightfile.txt'
		scoreTree_file = 'scoring_tree.tre'
		self.weight_tree_builder.buildScoreTree(refMSA_file, scoreTree_file)
		ordered_weights = self.BMweights(scoreTree_file, weightfile)

		shutil.copy(weightfile, BootDir)
		shutil.copy(refaln_file, BootDir)
		shutil.copy(prealn_file, BootDir)
		os.chdir(BootDir)		
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn, refMSA_file, n, numprocesses)
		
		print "Scoring, BMweights"
		BMscores=self.scorer.scoreMSA_Weighted(refMSA_file, n, numseq, alnlen, pflag, ordered_weights, weightfile, finalscore_fileBM)

		return(numseq, alnlen, BMscores)
		
		

class pdBootstrapper(Bootstrapper):
	''' Score with patristic distance only.'''
	
	def __init__(self, aligner, tree_builder, weight_tree_builder, scorer):
		'''initialization function'''
		self.weight_tree_builder = weight_tree_builder
		super(patristicBootstrapper, self).__init__(aligner, tree_builder, scorer)
		
		
	def runBootstrap(self, BootDir, prealn, refMSA_file, pflag, n, numprocesses, allscores_file): 
		
		# Create the scoring tree
		dist_matrix_file='dist_matrix.txt'
		scoreTree_file = 'scoring_tree.tre'
		
		self.weight_tree_builder.buildScoreTree(refMSA_file, scoreTree_file)	
		dist_matrix = self.weight_tree_builder.PDweights(scoreTree_file, dist_matrix_file, numseq)
		
		shutil.copy(dist_matrix_file, BootDir)
		shutil.copy(refaln_file, BootDir)
		shutil.copy(prealn_file, BootDir)
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn, refMSA_file, n, numprocesses)
		
		print "Scoring, PDweights"
		PDscores = self.scorer.scoreMSA_Patristic(refMSA_file, n, numseq, alnlen, pflag, dist_matrix, dist_matrix_file, finalscore_filePD)
		
		return(numseq, alnlen, PDscores)
				
		
		
		
		
class AllBootstrapper(Bootstrapper):
	'''Score with all algorithms.'''
	
	def __init__(self, aligner, tree_builder, weight_tree_builder, scorer):
		'''initialization function'''
		self.weight_tree_builder = weight_tree_builder
		super(AllBootstrapper, self).__init__(aligner, tree_builder, scorer)	
		
	def runBootstrap(self, BootDir, prealn, refMSA_file, pflag, n, numprocesses, finalscore_fileG, finalscore_fileBM, finalscore_filePD): 
		
		# Create the scoring tree
		dist_matrix_file = 'dist_matrix_file.txt'
		weightfile = 'weightfile.txt'
		scoreTree_file = 'scoring_tree.tre'
		buildScoreTree(refMSA_file, scoreTree_file)
		ordered_weights = self.BMweights(scoreTree_file, weightfile)
		dist_matrix = self.weight_tree_builder.PDweights(scoreTree_file, dist_matrix_file, numseq)
		
		shutil.copy(weightfile, BootDir)
		shutil.copy(dist_matrix_file, BootDir)
		shutil.copy(refaln_file, BootDir)
		shutil.copy(prealn_file, BootDir)
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn, refMSA_file, n, numprocesses)
		
		
		# Conduct the scoring
		print "scoring Guidance"
		Gscores= self.scorer.scoreMSA_Guidance(refMSA_file, n, numseq, alnlen, pflag, finalscore_fileG)
		print "scoring BranchManager"
		BMscores=self.scorer.scoreMSA_Weighted(refMSA_file, n, numseq, alnlen, pflag, ordered_weights, weightfile, finalscore_fileBM)
		print "scoring Patristic"
		PDscores = self.scorer.scoreMSA_Patristic(refMSA_file, n, numseq, alnlen, pflag, dist_matrix, dist_matrix_file, finalscore_filePD)
		
		return(numseq, alnlen, Gscores, BMscores, PDscores)
			
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
