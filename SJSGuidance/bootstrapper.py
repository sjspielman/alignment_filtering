import subprocess, sys, shutil, os
from Bio import SeqIO
from numpy import *

class Bootstrapper(object): 
	def __init__(self, aligner, tree_builder, scorer):
		self.aligner = aligner
		self.tree_builder = tree_builder
		self.scorer = scorer
	
	def bootstrap( self, prealn, refaln_file, n, numprocesses):
		
		# Obtain reference alignment, as well as number of sequences and length of alignment
		refaln_seq=[]
		infile=open(refaln_file, 'r')
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
			
		# Separate into PROCESSED trees for given alignment software
		print "Formatting trees"
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
		
	def runBootstrap(self, BootDir, prealn, refaln_file, n, numprocesses, finalscore_fileG, finalscore_fileBM, finalscore_filePD, finalscore_fileG_penal, finalscore_fileBM_penal, finalscore_filePD_penal, weightTree_file, bmweights_file, pdweights_file): 
		
		shutil.copy('BranchManager.jar', BootDir)
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn, refaln_file, n, numprocesses)
		
		# Create the scoring tree
		(dist_matrix, ordered_bmweights) = self.weight_tree_builder.buildScoreTree(refaln_file, weightTree_file, bmweights_file, pdweights_file, numseq)
		
		# Conduct the scoring
		print "scoring Guidance"
		(Gscores, Gscores_P)= self.scorer.scoreMSA_Guidance(refaln_file, n, numseq, alnlen, finalscore_fileG, finalscore_fileG_penal)
		print "scoring BranchManager"
		(BMscores, BMscores_P)=self.scorer.scoreMSA_Weighted(refaln_file, n, numseq, alnlen, ordered_bmweights, bmweights_file, finalscore_fileBM, finalscore_fileBM_penal)
		print "scoring Patristic"
		(PDscores, PDscores_P) = self.scorer.scoreMSA_Patristic(refaln_file, n, numseq, alnlen, dist_matrix, pdweights_file, finalscore_filePD, finalscore_filePD_penal)
		
		return(numseq, alnlen, Gscores, BMscores, PDscores, Gscores_P, BMscores_P, PDscores_P)