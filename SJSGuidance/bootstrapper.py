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

		## this chunk is IQ-Tree
		#self.tree_builder.buildBootTrees(n, refaln_seq, numseq, alnlen, final_treefile)
		#print "Constructing INFERENCE trees"
		#self.tree_builder.buildTrees(n, final_treefile)
		
		## this chunk is unique FastTree
		# Create strictly unique bootstrapped trees
		#print "Creating STRICTLY UNIQUE bootstrapped trees"
		#temp_treefile="tempBStre.tre"
		#self.tree_builder.buildUniqueTrees(n, refaln_seq, numseq, alnlen, final_treefile, temp_treefile, bootseq)
		#print "done uniquing"
			
		# Separate into PROCESSED trees for given alignment software
		print "Formatting trees" # uses my line-by-line reimplementation of the ruby script
		self.aligner.processTrees(n, final_treefile) 	
		
		print "Building bootstrap alignments"
		self.aligner.multiMakeAlignmentsGT(prealn, n, numprocesses)
		
		return(numseq, alnlen)
		
		
		
		
class AllBootstrapper(Bootstrapper):
	''' Process with all three algorithms and both normalizations'''
	def __init__(self, aligner, tree_builder, weight_tree_builder, scorer):
		'''initialization function'''
		self.weight_tree_builder = weight_tree_builder
		super(AllBootstrapper, self).__init__(aligner, tree_builder, scorer)	
		
	def runBootstrap(self, BootDir, prealn, refMSA_file, n, numprocesses, finalscore_fileG, finalscore_fileBM, finalscore_filePD, finalscore_fileG_penal, finalscore_fileBM_penal, finalscore_filePD_penal, scoreTree_file, weightfile, dist_matrix_file): 
		
		shutil.copy('BranchManager.jar', BootDir)
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn, refMSA_file, n, numprocesses)
		
		# Create the scoring tree
		(dist_matrix, ordered_weights) = self.weight_tree_builder.buildScoreTree(refMSA_file, scoreTree_file, weightfile, dist_matrix_file, numseq)
		
		# Conduct the scoring
		print "scoring Guidance"
		(Gscores, Gscores_P)= self.scorer.scoreMSA_Guidance(refMSA_file, n, numseq, alnlen, finalscore_fileG, finalscore_fileG_penal)
		print "scoring BranchManager"
		(BMscores, BMscores_P)=self.scorer.scoreMSA_Weighted(refMSA_file, n, numseq, alnlen, ordered_weights, weightfile, finalscore_fileBM, finalscore_fileBM_penal)
		print "scoring Patristic"
		(PDscores, PDscores_P) = self.scorer.scoreMSA_Patristic(refMSA_file, n, numseq, alnlen, dist_matrix, dist_matrix_file, finalscore_filePD, finalscore_filePD_penal)
		
		return(numseq, alnlen, Gscores, BMscores, PDscores, Gscores_P, BMscores_P, PDscores_P)
			
		

class patristicBootstrapper(Bootstrapper):
	''' Process with PDweights, original normalization only'''
	def __init__(self, aligner, tree_builder, weight_tree_builder, scorer):
		'''initialization function'''
		self.weight_tree_builder = weight_tree_builder
		super(patristicBootstrapper, self).__init__(aligner, tree_builder, scorer)
		
		
	def runBootstrap(self, BootDir, prealn, refMSA_file, n, numprocesses, allscores_file, scoreTree_file): 
		
		# Create the scoring tree
		matrixfile='dist_matrix.txt'
		dist_matrix = self.weight_tree_builder.buildScoreTree(refMSA_file, scoreTree_file, scoreTree_file, matrixfile)	
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn, refMSA_file, n, numprocesses)
		
		print "Scoring using the PATRISTIC"
		pscores = self.scorer.scoreMSA_patristic(refMSA_file, n, numseq, alnlen, dist_matrix, allscores_file)
		return(numseq, alnlen, pscores)
		



class guidanceBootstrapper(Bootstrapper):
	''' Process with Guidance, original normalization only'''
	def __init__(self, aligner, tree_builder, scorer):
		'''initialization function'''
		super(guidanceBootstrapper, self).__init__(aligner, tree_builder, scorer)

	def runBootstrap(self, BootDir, prealn, refMSA_file, n, numprocesses, allscores_file):
		os.chdir(BootDir)
		(numseq, alnlen) = self.bootstrap(prealn, refMSA_file, n, numprocesses)
		print "Scoring using the unweighted (original guidance) algorithm"
		gscores= self.scorer.scoreMSA_Guidance(refMSA_file, n, numseq, alnlen, allscores_file)
		return(numseq, alnlen, gscores)

class bmBootstrapper(Bootstrapper):
	''' Process with BMweights, original normalization only'''
	def __init__(self, aligner, tree_builder, weight_tree_builder, scorer):
		'''initialization function'''
		self.weight_tree_builder = weight_tree_builder
		super(bmBootstrapper, self).__init__(aligner, tree_builder, scorer)
		
	def runBootstrap(self, BootDir, prealn, refMSA_file, n, numprocesses, allscores_file, scoreTree_file, weightfile): 
		
		# Create the scoring tree
		ordered_weights = self.weight_tree_builder.buildScoreTree(refMSA_file, scoreTree_file, weightfile)
		shutil.copy(weightfile, BootDir)
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn, refMSA_file, n, numprocesses)
		
		print "Scoring using the weighted algorithm"
		gweightscores=self.scorer.scoreMSA_originalWeighted(refMSA_file, weightfile, n, numseq, alnlen, ordered_weights, allscores_file)

		return(numseq, alnlen, gweightscores)
		

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
