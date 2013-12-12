import subprocess, os, sys, re, csv, shutil
from dendropy import *
from dendropy import TaxonSet, Tree, TreeList
from Bio import SeqIO
from random import randint


class TreeBuilder:
	def __init__(self):
		'''initialization function'''
		return
		
		

class builderFastTree(TreeBuilder):

	def __init__(self, executable, options):
		'''"executable" is the path to the FastTree executable, and "options" are options to be handed to executable'''
		self.executable = executable
		self.options = options
	
	def buildBootTrees( self, num, refaln_seq, numseq, alnlen, outfile):
		''' Create bootstrap trees in FastTree. First, bootstrap the alignment 100x, then create a tree using those bootstrapped alignments.'''
		bootseq = 'refaln.BS'
		self.makeBootAlignments(num, refaln_seq, numseq, alnlen, bootseq)
		BuildTree=self.executable+' '+self.options+' -nosupport -n '+str(num)+' '+bootseq+' > '+outfile
		runtree=subprocess.call(str(BuildTree), shell='True')	
		return 0
		
	def makeBootAlignments(self, n, refaln_seq, numseq, alnlen, outfile):
		'''into a SINGLE FILE: makes n bootstrap alignment replicates from a given reference alignment sequence (refaln_seq). Doesn't involve alignment software at all.'''
		outhandle=open(outfile, 'w')
		for i in range(n):
			outhandle.write(' '+str(numseq)+' '+str(alnlen)+'\n')
			indices = []
			for a in range(alnlen):
				indices.append(randint(0,alnlen-1))	
			for s in range(numseq):
				newseq=''
				id=s+1
				for a in indices:
					newseq=newseq+refaln_seq[s][a]
				outhandle.write(str(id)+'        '+newseq+'\n')
		outhandle.close()
		
	def buildUniqueTrees(self, n, refaln_seq, numseq, alnlen, final_treefile, temp_treefile, bootseq):
		'''Constructs the bootstrap trees, but forcing that each tree is unique. The original Guidance did this it seemed, but only for mafft. 1. it's statistically insane and not actually how bootstrapping works, and 2. discrepancy=fail.'''
		baseTrees=[]
		redo=False
		numsametrees=0
		## Make a single bootstrap tree
		self.buildTrees(1, refaln_seq, numseq, alnlen, final_treefile, bootseq)
		atree = Tree(stream=open(final_treefile), schema="newick")
		baseTrees.append(atree)
		numtrees=1
		
		finalhandle=open(final_treefile, 'a')
		while numtrees<n:
			self.buildTrees(1, refaln_seq, numseq, alnlen, temp_treefile, bootseq)
			testTree = Tree(stream=open(temp_treefile), schema="newick")
			for baseTree in baseTrees:
				dist=baseTree.symmetric_difference(testTree)
				if dist == 0:
					numsametrees+=1
					redo=True
					break
				else:
					redo=False
			if redo==False:
				baseTrees.append(testTree)
				finalhandle.write(str(testTree)+'\n')
				numtrees+=1
				print numtrees, "so far"
		finalhandle.close()
		print "Found", str(numsametrees), "same trees."
		return 0			


class builderSemphy(TreeBuilder):
	'''Generate bootstrap trees with semphy, used by the original Guidance. Ours will use FastTree by default.'''
	def __init__(self, executable, options):
		self.executable = executable
		self.options = options
	
	def buildBootTrees(self, n, refaln_seq, numseq, alnlen, outfile):		
		''' Builds with semphy. DOESN'T need all those arguments, but FastTree does, and I'm trying to keep bootstrapper looking nice.'''
		BuildTree=self.executable+' '+self.options+' -s refaln.fasta -o out.out -T semphytree.tre -l log.out '
		runtree=subprocess.call(str(BuildTree), shell='True')	
		self.extractTrees(n, 'log.out', outfile)
		return 0
		
	def extractTrees(self, n, semphylog, outfile):
		treecounter=0
		inhandle=open(semphylog, 'r')
		lines=inhandle.readlines()
		inhandle.close()
		outhandle=open(outfile, 'w')
		numlines=len(lines)
		
		## Start collecting trees at line 9 (index 8).
		for a in range(8, numlines):
			tree=str(lines[a])
			outhandle.write(tree)
			treecounter+=1
			
			# Just in case. Likely unnecessary.		
			if treecounter==int(n):
				break
		outhandle.close()	
		return 0
			
