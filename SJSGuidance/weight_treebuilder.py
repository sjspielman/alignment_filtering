import subprocess, os, sys, re, csv, shutil
from dendropy import *
from dendropy import TaxonSet, Tree, TreeList
from Bio import SeqIO, AlignIO
from numpy import *

class WeightTreeBuilder:
	def __init__(self):
		'''initialization function'''
		return

	def killPolytomyDendro(self, infile):
		'''actually removes polytomies, except branchmanager is a beotch.'''
		rawtree = Tree(stream=open(infile), schema="newick")
		rawtree.resolve_polytomies(update_splits=True)
		newtree=str(rawtree)
		newtree.replace('[&U] ','')
		newtree.replace('[&R] ','')
		newtree2 = re.sub('\d+\.*\d*e-\d+', '0', newtree)
		out_handle = open(infile, 'w')
		out_handle.write(newtree2+';')
		out_handle.close()
		return 0
		
	
	def findWeights (self, treefile, weightfile):
		'''uses BranchManager to get phylogenetic weights and processes output for later use in scoring function'''
		callBM = 'java -cp BranchManager.jar BM '+treefile+' > '+weightfile
		runBM = subprocess.call(str(callBM), shell='True')
		
		# Process weights for subsequent use in scoring function
		handle=open(weightfile,'r')
		weights=handle.readlines()
		handle.close()
		weights_dict={}
		for line in weights:
			find = re.search('(\d+)\t(0\.\d+)', line)
			if find:
				taxon = int(find.group(1))
				weight = float(find.group(2))
				weights_dict[taxon]=weight
		ordered_weights=[]
		handle = open(weightfile, 'w')
		for key in sorted(weights_dict.iterkeys()):
			handle.write("%s\n" % (weights_dict[key]))
			ordered_weights.append(weights_dict[key])
		handle.close()
		
		return ordered_weights		
		
	def getPatristic(self, treefile, matrixfile, numseq):
		tree=Tree.get_from_path(treefile, 'newick')
		patmat=treecalc.PatristicDistanceMatrix(tree)
		pat_dict={}
		for i, t1 in enumerate(tree.taxon_set):
			for t2 in tree.taxon_set[i:]:
				# Smaller taxon name comes first in key
				if int(t1.label) < int(t2.label):
				   key=str(t1.label)+'_'+str(t2.label)
				else:
					key=str(t2.label)+'_'+str(t1.label)
				pat_dict[key]=patmat(t1,t2)
		
		## Normalizes as they get added into the matrix. 
		maxdist = pat_dict[max(pat_dict, key=pat_dict.get)]
		maxdist=1
		dist_matrix = zeros((numseq, numseq))
		for i in range(numseq):
			for x in range(numseq):
				if i==x:
					dist_matrix[i][x] = (0)
				elif i < x:
					dist_matrix[i][x] = pat_dict[str(i+1)+'_'+str(x+1)] / maxdist
				else:
					dist_matrix[i][x] = pat_dict[str(x+1)+'_'+str(i+1)] / maxdist
		savetxt(matrixfile, dist_matrix, delimiter=' ', fmt='%.5f')
		return dist_matrix
	
	
	
class bothRAxML(WeightTreeBuilder):	
	'''Does both branchmananger weights and patristic distance matrix.'''	
	def __init__(self, executable, options):
		self.executable = executable
		self.options = options
		
	def buildScoreTree(self, alnfile, treefile, weightfile, matrixfile, numseq):
		print "Building Scoring Tree with RAxML using amino acid data"
		
		## Convert alignment to phylip for fasta input
		aln = AlignIO.read(alnfile, 'fasta')
		AlignIO.write(aln, 'temp.phy', 'phylip')
		
		BuildTree=self.executable+' '+self.options+' -s temp.phy -n out'
		subprocess.call(BuildTree, shell=True)
		shutil.move('RAxML_bestTree.out', treefile)
		###subprocess.call('rm RAxML*', shell=True)
		
		dist_matrix= self.getPatristic(treefile, matrixfile, numseq)
		ordered_weights = self.findWeights(treefile, weightfile)
		
		return (dist_matrix, ordered_weights)
	
		
class bothFastTree(WeightTreeBuilder):	
	'''Does both branchmananger weights and patristic distance matrix.'''	
	def __init__(self, executable, options):
		self.executable = executable
		self.options = options

	def buildScoreTree(self, alnfile, treefile, weightfile, matrixfile, numseq):
		print "Building Scoring Tree with FastTree (-slow)"
		BuildTree=self.executable+' '+self.options+' -nosupport '+alnfile+' > '+treefile
		subprocess.call(BuildTree, shell=True)
		
		dist_matrix= self.getPatristic(treefile, matrixfile, numseq)
		ordered_weights = self.findWeights(treefile, weightfile)
		
		return (dist_matrix, ordered_weights)
		
				
		
class bothIQTree(WeightTreeBuilder):
	'''Does both branchmananger weights and patristic distance matrix.'''	
	def __init__(self, executable, options):
		self.executable = executable
		self.options = options

	def buildScoreTree(self, alnfile, treefile, weightfile, matrixfile, numseq):
		print "Building Scoring Tree with IQ-Tree. (default model=WAG)"
		BuildTree=self.executable+' '+self.options+' -s '+alnfile
		subprocess.call(BuildTree, shell=True)
		shutil.move(alnfile+'.treefile', treefile)
		
		dist_matrix= self.getPatristic(treefile, matrixfile, numseq)

		ordered_weights = self.findWeights(treefile, weightfile)
		return (dist_matrix, ordered_weights)


				

class patIQTree(WeightTreeBuilder):
	def __init__(self, executable, options):
		self.executable = executable
		self.options = options
	
	def buildScoreTree(self, alnfile, treefile, matrixfile, numseq):
		print "Building Scoring Tree with IQ-Tree. (default model=WAG)"
		BuildTree=self.executable+' '+self.options+' -s '+alnfile
		subprocess.call(BuildTree, shell=True)
		shutil.move(alnfile+'.treefile', treefile)
		pat_dict = self.getPatristic(treefile, matrixfile, numseq)
		return pat_dict

class patFastTree(WeightTreeBuilder):
	def __init__(self, executable, options):
		self.executable = executable
		self.options = options
	
	def buildScoreTree(self, alnfile, treefile, matrixfile):
		print "Building Scoring Tree"
		BuildTree=self.executable+' '+self.options+' '+alnfile+' > '+treefile
		runtree=subprocess.call(str(BuildTree), shell='True')	
		print "Computing Distances"
		dist_matrix = self.getPatristic(treefile, matrixfile)
		return dist_matrix
		

class weightTreeFasttree(WeightTreeBuilder):
	def __init__(self, executable, options):
		self.executable = executable
		self.options = options
		
	def buildScoreTree(self, alnfile, treefile, weightfile):
		print "Building Scoring Tree"
		BuildTree=self.executable+' '+self.options+' '+alnfile+' > '+treefile
		runtree=subprocess.call(str(BuildTree), shell='True')
		print "Processing Polytomies (required by BranchManager)"
		self.killPolytomyDendro(treefile)	
		print "Computing Taxon Weights"
		ordered_weights = self.findWeights(treefile, weightfile)
		
		return ordered_weights

class weightIQTree(WeightTreeBuilder):
	def __init__(self, executable, options):
		self.executable = executable
		self.options = options
	
	def buildScoreTree(self, alnfile, treefile, weightfile):
		print "Building Scoring Tree with IQ-Tree. (default model=WAG)"
		BuildTree=self.executable+' '+self.options+' -s '+alnfile
		subprocess.call(BuildTree, shell=True)
		shutil.move(alnfile+'.treefile', treefile)
		
		print "Processing Polytomies (required by BranchManager)"
		self.killPolytomyDendro(treefile)	
		print "Computing Taxon Weights"
		ordered_weights = self.findWeights(treefile, weightfile)
		
		# Remove iqtree files
		rmiq='rm '+alnfile+'.*'
		subprocess.call(rmiq, shell=True)
		
		return ordered_weights
		
		
class weightTreeRAXML(WeightTreeBuilder):
	def __init__(self, executable, options):
		self.executable = executable
		self.options = options
		
	def buildScoreTree(self, alnfile, treefile, weightfile):
		# Convert alnfile into phylip for the RAxML builder
		parsed=AlignIO.read(alnfile, 'fasta')
		AlignIO.write(parsed, 'refaln.phy', 'phylip')
		
		print "Building Scoring Tree (this may be quite slow with RAxML!)"
		BuildTree=self.executable+' '+self.options+'-d -p 6957 -s refaln.phy -n finaltree'
		print BuildTree
		runtree=subprocess.call(str(BuildTree), shell='True')
		
		# Move RAxML tree to the treefile
		command = "mv RAxML_result.finaltree "+treefile
		subprocess.call(command, shell=True)
		
		print "Processing Polytomies (required by BranchManager)"
		self.killPolytomyDendro(treefile)	
		print "Computing Taxon Weights"
		ordered_weights = self.findWeights(treefile, weightfile)
		
		print "Removing RAxML vomit"
		command = 'rm RAxML_*'
		subprocess.call(command, shell=True)
		
		return ordered_weights
		
		
		
		
		
		
		
		





