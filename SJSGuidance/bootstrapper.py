import subprocess, sys, shutil, os
from Bio import SeqIO
from numpy import *

class Bootstrapper(object): 
	def __init__(self, aligner, tree_builder):
		'''initialization function'''
		self.aligner = aligner
		self.tree_builder = tree_builder
	
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
	
	def __init__(self, aligner, tree_builder):
		'''initialization function'''
		super(guidanceBootstrapper, self).__init__(aligner, tree_builder)

	def runBootstrap(self, BootDir, prealn, refMSA_file, pflag, n, numprocesses):
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn, refMSA_file, n, numprocesses)

		return(numseq, alnlen)


class bmBootstrapper(Bootstrapper):
	''' Score with BranchManager weights only.'''
	
	def __init__(self, aligner, tree_builder, weight_tree_builder):
		'''initialization function'''
		self.weight_tree_builder = weight_tree_builder
		super(bmBootstrapper, self).__init__(aligner, tree_builder)
		
	def runBootstrap(self, BootDir, prealn, refMSA_file, pflag, n, numprocesses): 
		
		# Create the scoring tree
		weightfile='weightfile.txt'
		scoreTree_file = 'scoring_tree.tre'
		self.weight_tree_builder.buildScoreTree(refMSA_file, scoreTree_file)
		ordered_weights = self.weight_tree_builder.calcBMweights(scoreTree_file, weightfile)

		shutil.copy(weightfile, BootDir)
		shutil.copy(refaln_file, BootDir)
		shutil.copy(prealn_file, BootDir)
		os.chdir(BootDir)		
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn, refMSA_file, n, numprocesses)

		return(numseq, alnlen)
		
		

class pdBootstrapper(Bootstrapper):
	''' Score with patristic distance only.'''
	
	def __init__(self, aligner, tree_builder, weight_tree_builder):
		'''initialization function'''
		self.weight_tree_builder = weight_tree_builder
		super(patristicBootstrapper, self).__init__(aligner, tree_builder)
		
		
	def runBootstrap(self, BootDir, prealn, refMSA_file, pflag, n, numprocesses): 
		
		# Create the scoring tree
		dist_matrix_file='dist_matrix.txt'
		scoreTree_file = 'scoring_tree.tre'
		
		self.weight_tree_builder.buildScoreTree(refMSA_file, scoreTree_file)	
		dist_matrix = self.weight_tree_builder.calcPDweights(scoreTree_file, dist_matrix_file, numseq)
		
		shutil.copy(dist_matrix_file, BootDir)
		shutil.copy(refaln_file, BootDir)
		shutil.copy(prealn_file, BootDir)
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn, refMSA_file, n, numprocesses)

		return(numseq, alnlen)
				
		
		
		
		
class AllBootstrapper(Bootstrapper):
	'''Score with all algorithms.'''
	
	def __init__(self, aligner, tree_builder, weight_tree_builder):
		'''initialization function'''
		self.weight_tree_builder = weight_tree_builder
		super(AllBootstrapper, self).__init__(aligner, tree_builder)	
		
	def runBootstrap(self, BootDir, prealn, refMSA_file, pflag, n, numprocesses): 
		
		# Create the scoring tree
		dist_matrix_file = 'dist_matrix_file.txt'
		weightfile = 'weightfile.txt'
		scoreTree_file = 'scoring_tree.tre'
		buildScoreTree(refMSA_file, scoreTree_file)
		ordered_weights = self.weight_tree_builder.calcBMweights(scoreTree_file, weightfile)
		dist_matrix = self.weight_tree_builder.calcPDweights(scoreTree_file, dist_matrix_file, numseq)
		
		shutil.copy(weightfile, BootDir)
		shutil.copy(dist_matrix_file, BootDir)
		shutil.copy(refaln_file, BootDir)
		shutil.copy(prealn_file, BootDir)
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn, refMSA_file, n, numprocesses)
	
		return(numseq, alnlen)
			