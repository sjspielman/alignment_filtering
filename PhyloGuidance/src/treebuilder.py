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
	
	def makeBootAlignment(self, refaln_seq, numseq, alnlen, outfile):
		'''into a SINGLE FILE: makes n bootstrap alignment replicates from a given reference alignment sequence (refaln_seq). Doesn't involve alignment software at all.'''
		outhandle=open(outfile, 'w')
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
	

	def buildBootTree( self, refaln_seq, numseq, alnlen, outfile):
		bootseq = 'refaln.BS'
		self.makeBootAlignment(refaln_seq, numseq, alnlen, bootseq)
		BuildTree=self.executable+' '+self.options+' -nosupport '+bootseq+' > '+outfile
		runtree=subprocess.call(str(BuildTree), shell='True')	
		
		# Double-check that FastTree worked. (It has been throwing some floating point exceptions.) If not, make a new bootstrap alignment and try again.
		numReps = 0
		while os.path.getsize(outfile) <= 0:
			self.makeBootAlignments(num, refaln_seq, numseq, alnlen, bootseq)
			runtree=subprocess.call(str(BuildTree), shell='True')	
			assert(numReps<10), "Serious FastTree problem."
			numReps+=1


	def buildBootTreesNoReps(self, num, refaln_seq, numseq, alnlen, outfile):
		''' Construct bootstrap trees, but if found a repeat just remember which ones are repeats and save that many alns from it later '''
		
		saveTrees = []
		numSaveTrees = []
		
		## Final file
		finalTrees = open(outfile, 'w')
	
		self.buildBootTree(refaln_seq, numseq, alnlen, 'temp.tre')
		
		# Save tree
		tree = Tree(stream=open('temp.tre'), schema="newick")
		finalTrees.write(str(tree)+';\n')
		saveTrees.append(tree)
		numSaveTrees.append(1)
		
		for i in range(1,num):
			self.buildBootTree(refaln_seq, numseq, alnlen, 'temp.tre')
			
			# Compare it to all current trees to see if we already have it
			testTree = Tree(stream=open('temp.tre'), schema="newick")
			newTree = True
			for x in range(len(saveTrees)):
				dist = saveTrees[x].symmetric_difference(testTree)
				
				# Have the tree already. Increment 1 to that tree's index
				if ( dist - 1. < 1e-10 ):
					numSaveTrees[x]+=1
					newTree = False
					break
			
			# Need to save the tree
			if newTree:
				numSaveTrees.append(1)
				saveTrees.append(testTree)
				finalTrees.write(str(testTree)+';\n')
		
		finalTrees.close()
		
		# Double check that the correct number of trees (num bootstraps requested) have been accounted for
		sum = 0
		for i in numSaveTrees:
			sum += i
		assert (sum == num), "The correct number of trees have not been built. This is a problem. Email stephanie.spielman@gmail.com "
		
		return numSaveTrees
		
		
			

#####################################################################################################################
#####################################################################################################################
####### FOLLOWING CLASSES ARE NOT USED AND ARE PROBABLY SOMEWHAT DEPRECATED. YOU CAN REVIVE THEM AT YOUR OWN RISK!!!!
####################################################################################################################

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
			