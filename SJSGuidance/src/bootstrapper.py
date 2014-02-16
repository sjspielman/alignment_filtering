import subprocess, sys, shutil, os
from Bio import SeqIO
from numpy import *

class Bootstrapper(object): 
	def __init__(self, aligner, tree_builder, scorer):
		self.aligner = aligner
		self.tree_builder = tree_builder
		self.scorer = scorer
	
	def bootstrap( self, prealn_file, refaln_file, n, numprocesses):
		
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
		self.aligner.multiMakeAlignmentsGT(prealn_file, n, numprocesses)
		
		return(numseq, alnlen)
		
		
		
		
class AllBootstrapper(Bootstrapper):
	''' Process with all three algorithms and both normalizations'''
	def __init__(self, aligner, tree_builder, weight_tree_builder, scorer):
		'''initialization function'''
		self.weight_tree_builder = weight_tree_builder
		super(AllBootstrapper, self).__init__(aligner, tree_builder, scorer)	
		
	def runBootstrap(self, BootDir, prealn_file, refaln_file, n, numprocesses, weightTree_file, bmweights_file, pdweights_file): 
		
		
		shutil.copy('src/BranchManager.jar', BootDir)
		shutil.copy(refaln_file, BootDir)
		shutil.copy(prealn_file, BootDir)
		os.chdir(BootDir)
		
		# Call bootstrapper
		(numseq, alnlen) = self.bootstrap(prealn_file, refaln_file, n, numprocesses)
		
		# Create the scoring tree
		(dist_matrix, ordered_bmweights) = self.weight_tree_builder.buildScoreTree(refaln_file, weightTree_file, bmweights_file, pdweights_file, numseq)
		
		
		## Final score files names
		g="scores_Guidance.txt"
		bm="scores_BMweights.txt"
		pd="scores_PDweights.txt"
		gP="scores_GuidanceP.txt"
		bmP="scores_BMweightsP.txt"
		pdP="scores_PDweightsP.txt"
		
		
		# Conduct the scoring
		print "scoring Guidance"
		(gscores, gscores_p)= self.scorer.scoreMSA_Guidance(refaln_file, n, numseq, alnlen, g, gP)
		print "scoring BranchManager"
		(bmscores, bmscores_p)=self.scorer.scoreMSA_BMweights(refaln_file, n, numseq, alnlen, ordered_bmweights, bmweights_file, bm, bmP)
		print "scoring Patristic"
		(pdscores, pdscores_p) = self.scorer.scoreMSA_PDweights(refaln_file, n, numseq, alnlen, dist_matrix, pdweights_file, pd, pdP)
		
		# Place scores into dictionary. Useful for naming final files.
		alg_scores={'Guidance':gscores, 'BMweights':bmscores, 'PDweights':pdscores, 'GuidanceP':gscores_p, 'BMweightsP':bmscores_p, 'PDweightsP':pdscores_p}

		return(numseq, alnlen, alg_scores)
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		