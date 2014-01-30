import subprocess
import re
import shutil
from dendropy import Tree, treecalc
from Bio import AlignIO


class WeightTreeBuilder:
	def __init__(self):
		'''initialization function'''
		return

	def rmPolytomy(self, infile):
		'''Changes polytomies to branch length=0 using DendroPy. Useful for BranchManager, but we should be consistent and use for all weighted scoring.'''
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
		
	
	def calcBMweights(self, treefile, weightfile):
		'''uses BranchManager to get phylogenetic weights and processes output for later use in scoring function'''
		
		self.rmPolytomy(treefile)
		
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
		
	def calcPDweights(self, treefile, matrixfile, numseq):
		
		self.rmPolytomy(treefile)
		
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
		

	
class weightTreeRAxML(WeightTreeBuilder):	
	'''Build phylogeny from which to calculate phylogenetic weights using RAxML'''	
	def __init__(self, executable, options):
		self.executable = executable
		self.options = options
		
	def buildScoreTree(self, alnfile, treefile):
		print "Building Weighting Tree with RAxML using amino acid data"
		
		## Convert alignment to phylip for fasta input
		aln = AlignIO.read(alnfile, 'fasta')
		AlignIO.write(aln, 'temp.phy', 'phylip')
		
		BuildTree=self.executable+' '+self.options+' -s temp.phy -n out'
		subprocess.call(BuildTree, shell=True)
		shutil.move('RAxML_result.out', treefile)
		subprocess.call('rm RAxML*', shell=True) #Removes the extra files RAxML produces when constructing phylogenies
	
		return 0
	
	
	
class weightTreeFastTree(WeightTreeBuilder):	
	'''Build phylogeny from which to calculate phylogenetic weights using FastTree'''	
	def __init__(self, executable, options):
		self.executable = executable
		self.options = options

	def buildScoreTree(self, alnfile, treefile):
		print "Building Weighting Tree with FastTree"
		BuildTree=self.executable+' '+self.options+' -nosupport '+alnfile+' > '+treefile
		subprocess.call(BuildTree, shell=True)
		
		return 0

