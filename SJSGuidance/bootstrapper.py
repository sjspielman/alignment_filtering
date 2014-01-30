import subprocess, sys, shutil, os
from Bio import SeqIO
from numpy import *

class Bootstrapper(object): 
	def __init__(self, aligner, tree_builder):
		'''initialization function'''
		self.aligner = aligner
		self.tree_builder = tree_builder
	
	def bootstrap( self, BootDir, prealn_file, refMSA_file, n, numprocesses):
		
		os.chdir(BootDir)
		
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
		self.aligner.multiMakeAlignmentsGT(prealn_file, n, numprocesses)
		
		return(numseq, alnlen)
		


class guidanceBootstrapper(Bootstrapper):
	''' Score with original guidance only.'''
	
	def __init__(self, aligner, tree_builder):
		'''initialization function'''
		super(guidanceBootstrapper, self).__init__(aligner, tree_builder)

	def runBootstrap(self, BootDir, prealn_file, refMSA_file, pflag, n, numprocesses):
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn_file, refMSA_file, n, numprocesses)

		return(numseq, alnlen)


class bmBootstrapper(Bootstrapper):
	''' Score with BranchManager weights only.'''
	
	def __init__(self, aligner, tree_builder, weight_tree_builder):
		'''initialization function'''
		self.weight_tree_builder = weight_tree_builder
		super(bmBootstrapper, self).__init__(aligner, tree_builder)
		
	def runBootstrap(self, BootDir, prealn_file, refMSA_file, pflag, n, numprocesses): 
		
		# Create the scoring tree. Remove polytomies from scoring tree as BranchManager requires fully bifurcating trees
		weightfile='weightfile.txt'
		weightTree_file = 'weighting_tree.tre'
		self.weight_tree_builder.buildScoreTree(refMSA_file, weightTree_file)
		self.weight_tree_builder.rmPolytomy(weightTree_file)
		ordered_weights = self.weight_tree_builder.calcBMweights(weightTree_file, weightfile)

		shutil.copy(weightfile, BootDir)
		shutil.copy(refaln_file, BootDir)
		shutil.copy(prealn_file_file, BootDir)
		os.chdir(BootDir)		
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn_file, refMSA_file, n, numprocesses)

		return(numseq, alnlen)
		
		

class pdBootstrapper(Bootstrapper):
	''' Score with patristic distance only.'''
	
	def __init__(self, aligner, tree_builder, weight_tree_builder):
		'''initialization function'''
		self.weight_tree_builder = weight_tree_builder
		super(PDBootstrapper, self).__init__(aligner, tree_builder)
		
		
	def runBootstrap(self, BootDir, prealn_file, refMSA_file, pflag, n, numprocesses): 
		
		# Create the scoring tree
		dist_matrix_file='dist_matrix.txt'
		weightTree_file = 'weighting_tree.tre'
		
		self.weight_tree_builder.buildScoreTree(refMSA_file, weightTree_file)	
		dist_matrix = self.weight_tree_builder.calcPDweights(weightTree_file, dist_matrix_file, numseq)
		
		shutil.copy(dist_matrix_file, BootDir)
		shutil.copy(refaln_file, BootDir)
		shutil.copy(prealn_file_file, BootDir)
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn_file, refMSA_file, n, numprocesses)

		return(numseq, alnlen)
				
		
		
		
		
class WeightedBootstrapper(Bootstrapper):
	'''Bootstrap for both BMweights and PDweights.'''
	
	def __init__(self, aligner, tree_builder, weight_tree_builder):
		'''initialization function'''
		self.weight_tree_builder = weight_tree_builder
		super(WeightedBootstrapper, self).__init__(aligner, tree_builder)	
		
	def runBootstrap(self, BootDir, prealn_file, refMSA_file, pflag, n, numprocesses): 
		
		# Create the scoring tree. Remove polytomies as BranchManager requires fully bifurcating. As both BMweights and PDweights here, should be consistent and use the same non-polytomy tree for both
		dist_matrix_file = 'dist_matrix_file.txt'
		weightfile = 'weightfile.txt'
		weightTree_file = 'weighting_tree.tre'
		
		self.weight_tree_builder.buildScoreTree(refMSA_file, weightTree_file)
		self.weight_tree_builder.rmPolytomy(weightTree_file)
		ordered_weights = self.weight_tree_builder.calcBMweights(weightTree_file, weightfile)
		dist_matrix = self.weight_tree_builder.calcPDweights(weightTree_file, dist_matrix_file, numseq)
		
		shutil.copy(weightfile, BootDir)
		shutil.copy(dist_matrix_file, BootDir)
		shutil.copy(refaln_file, BootDir)
		shutil.copy(prealn_file_file, BootDir)
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn_file, refMSA_file, n, numprocesses)
	
		return(numseq, alnlen)
			