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
	

	def buildBootTrees( self, num, refaln_seq, numseq, alnlen, outfile):
		bootseq = 'refaln.BS'
		self.makeBootAlignments(num, refaln_seq, numseq, alnlen, bootseq)
		BuildTree=self.executable+' '+self.options+' -nosupport -n '+str(num)+' '+bootseq+' > '+outfile
		runtree=subprocess.call(str(BuildTree), shell='True')	
		return 0


	## This function is currently not used, and we do not advocate its use. It was originally implemented to fully re-implement Guidance, but we do not feel that forcing bootstrap replicates to be unique accurately reflects the spirit of the bootstrap.
	def buildUniqueTrees(self, n, refaln_seq, numseq, alnlen, final_treefile, temp_treefile, bootseq):
		'''Constructs the bootstrap trees, but forcing that each tree is unique.'''
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

###########################################################################################################
###########################################################################################################
####### FOLLOWING CLASSES ARE NOT USED, BUT YOU ARE WELCOME TO USE THEM IF YOU WANT TO! sjs, 2/16/14 ######
###########################################################################################################

## Note that if you want to use this, you'll need to set the treebuilder class in main.py as follows - (but double check your executable path)
###     tmod=builderSemphy("../semphy/semphy", " -a 20 --jtt -H -J -v 5 --BPrepeats="+str(n)+" ")
class builderSemphy(TreeBuilder):
	def __init__(self, executable, options):
		self.executable = executable
		self.options = options
	
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
		
	def buildTrees(self, n, refaln_seq, numseq, alnlen, final_treefile, bootseq):		
		''' Builds with semphy. DOESN'T need all those arguments, but FastTree does, and I'm trying to keep bootstrapper looking nice.'''
		BuildTree=self.executable+' '+self.options+' -s refaln.fasta -o out.out -T semphytree.tre -l log.out '
		runtree=subprocess.call(str(BuildTree), shell='True')	
		self.extractTrees(n, 'log.out', final_treefile)
		return 0


class builderIQTree(TreeBuilder):
	def __init__(self, executable, options):
		self.executable = executable
		self.options = options	

	def buildTrees(self, n, final_treefile):		
		''' Builds 100 tree inferences [not bootstrapped] with iqtree. DOESN'T need all those arguments, but FastTree does, and I'm trying to keep bootstrapper looking nice.'''
		BuildTree=self.executable+' '+self.options+' -s refaln.fasta -n '+str(n)+' -wt'
		runtree=subprocess.call(str(BuildTree), shell='True')	
		shutil.move('refaln.fasta.treels', final_treefile)
		return 0
	
	
	def buildBootTrees(self, n, refaln_seq, numseq, alnlen, final_treefile):
		'''Bootstraps the sequence (each into its own file!!) and performs a single inference on each bootstrap.'''
		
		bsaln_file='bsaln.fasta'
		treefile='bsaln.fasta.treefile'
		BStreefile=open(final_treefile, 'a')
		
		for i in range(n):
			indices=[]
			outhandle=open(bsaln_file, 'w')
			for a in range(alnlen):
				indices.append(randint(0,alnlen-1))	
			for s in range(numseq):
				newseq=''
				id=s+1
				for a in indices:
					newseq=newseq+refaln_seq[s][a]
				outhandle.write('>'+str(id)+'\n'+newseq+'\n')
			outhandle.close()
			
			BuildTree=self.executable+' '+self.options+' -s '+bsaln_file+' -n 1' 
			runtree=subprocess.call(str(BuildTree), shell='True')	
			treeh = open(treefile, 'r')
			tree=treeh.read()
			treeh.close()
			
			BStreefile.write(tree+'\n')

		BStreefile.close()
		
		command='rm bsaln*'
		subprocess.call(command, shell=True)
		return 0
			