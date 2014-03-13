import subprocess, sys, shutil, os
from Bio import SeqIO
from numpy import *
import multiprocessing

class Bootstrapper(object): 
	def __init__(self, **kwargs):
		self.aligner      = kwargs.get("aligner")
		self.tree_builder = kwargs.get("tree_builder")
		self.scorer       = kwargs.get("scorer")
		
		self.n            = kwargs.get("bootstraps")
		self.prealn_file  = kwargs.get("prealn_file")
		self.refaln_file  = kwargs.get("refaln_file")
		self.final_treefile = "BStrees.tre"
		self.BootDir      = kwargs.get("BootDir", "BootDir/")
		
		# Don't overload in case 
		availableCPU = multiprocessing.cpu_count() 
		if availableCPU > 1:
			availableCPU -= 1
			
		self.numprocesses = kwargs.get("threads", availableCPU)
		print "cpu", multiprocessing.cpu_count()
		
		############# input assertions ##########
		assert(self.aligner is not None), "No aligner was passed to Bootstrapper."
		assert(self.tree_builder is not None), "No tree builder was passed to Bootstrapper."
		assert(self.scorer is not None), "No scorer was passed to Bootstrapper."
		assert(self.n is not None), "Number of bootstraps was not specified."
		assert(self.prealn_file is not None), "Raw sequence file was not specified."
		assert(self.refaln_file is not None), "Reference alignment file was not specified."
		assert(self.numprocesses is not None), "You appear to have no CPUs on your system."
		#########################################

		# Parse reference alignment for internal use
		self.refaln_seq=[]
		infile=open(self.refaln_file, 'r')
		parsed = list(SeqIO.parse(infile, 'fasta'))
		infile.close()
		for record in parsed:
			self.refaln_seq.append(str(record.seq))
		self.numseq = len(self.refaln_seq)
		self.alnlen = len(self.refaln_seq[0])	
		
		#Final file name information
		
		
	def bootstrap(self):
		''' Create bootstrapped trees and then from those create perturbed alignments.'''
		
		# Create the bootstrapped trees
		print "Constructing bootstrap trees"
		self.tree_builder.buildBootTrees(self.n, self.refaln_seq, self.numseq, self.alnlen, self.final_treefile)
			
		# Separate into PROCESSED trees for given alignment software
		print "Formatting trees"
		self.aligner.processTrees(self.n, self.final_treefile) 	
		
		print "Building bootstrap alignments"
		self.aligner.multiMakeAlignmentsGT(self.prealn_file, self.n, self.numprocesses)
				
		
class AllBootstrapper(Bootstrapper):
	''' Process with all three algorithms and both normalizations'''
	def __init__(self, **kwargs):
		super(AllBootstrapper, self).__init__(**kwargs)	
		self.weight_tree_builder = kwargs.get("weight_tree_builder")
		assert(self.weight_tree_builder is not None), "No weight tree builder was passed to Bootstrapper."

		# Names for these files don't really matter.
		self.weightTree_file = kwargs.get("weightTree_file", "scoringtree.tre")
		self.bmweights_file  = kwargs.get("bmweights_file", "bmweights.txt")
		self.pdweights_file  = kwargs.get("pdweights_file", "pd_matrix.txt")
		
	def runBootstrap(self): 
			
		shutil.copy('src/BranchManager.jar', self.BootDir)
		shutil.copy(self.refaln_file, self.BootDir)
		shutil.copy(self.prealn_file, self.BootDir)
		os.chdir(self.BootDir)
		
		# Call bootstrapper
		self.bootstrap()
		
		# Create the scoring tree
		(dist_matrix, ordered_bmweights) = self.weight_tree_builder.buildScoreTree(self.refaln_file, self.weightTree_file, self.bmweights_file, self.pdweights_file, self.numseq)
		
		
		## Final score files names
		g="scores_Guidance.txt"
		bm="scores_BMweights.txt"
		pd="scores_PDweights.txt"
		gP="scores_GuidanceP.txt"
		bmP="scores_BMweightsP.txt"
		pdP="scores_PDweightsP.txt"
		
		
		# Conduct the scoring
		print "scoring Guidance"
		(gscores, gscores_p)= self.scorer.scoreMSA_Guidance(self.refaln_file, self.n, self.numseq, self.alnlen, g, gP)
		print "scoring BranchManager"
		(bmscores, bmscores_p)=self.scorer.scoreMSA_BMweights(self.refaln_file, self.n, self.numseq, self.alnlen, ordered_bmweights, self.bmweights_file, bm, bmP)
		print "scoring Patristic"
		(pdscores, pdscores_p) = self.scorer.scoreMSA_PDweights(self.refaln_file, self.n, self.numseq, self.alnlen, dist_matrix, self.pdweights_file, pd, pdP)
		
		# Place scores into dictionary. Useful for naming final files.
		alg_scores={'Guidance':gscores, 'BMweights':bmscores, 'PDweights':pdscores, 'GuidanceP':gscores_p, 'BMweightsP':bmscores_p, 'PDweightsP':pdscores_p}

		return(self.numseq, self.alnlen, alg_scores)
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		