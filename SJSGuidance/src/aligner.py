import subprocess, sys, re, os, multiprocessing
from multiprocessing import pool
from dendropy import *
from Bio import AlignIO
    
class Aligner:
	def __init__(self):
		return

	def makeAlignment( self, prealn_file, alnfile ):
		'''Makes an alignment, nothing fancy.'''
		print "Aligner base class. Function makeAlignment() is not implemented."

	def makeAlignmentGT( self, treefile, prealn_file, alnfile):
		'''Makes an alignment using a provided guide tree.'''
		print "Aligner base class. Function makeAlignmentGT() is not implemented."
	
	def multiMakeAlignmentsGT(self, prealn_file, n, numprocesses):
		'''Makes n bootstrap alignments. Note that this function is rather clunky because multiprocessing.pool() does not work when inside a class in python. The code here is a decent workaround of this unfortunate issue.'''
	
		if numprocesses=='':
			numprocesses=multiprocessing.cpu_count()	
		pool=multiprocessing.Pool(numprocesses)
		jobs=[]
		
		if numprocesses>=n:
			for i in range(n):
				treefile='tree'+str(i)+'.txt'
				alnfile='bootaln'+str(i)+'.fasta'
				p=multiprocessing.Process(target=self.makeAlignmentGT, args=( treefile, prealn_file, alnfile ))
				jobs.append(p)
				p.start()
			for p in jobs:
				p.join()
		
		elif numprocesses<n:
			nruns=0
			while (nruns < n):
				if (n-nruns)>=numprocesses:
					jobs=[]
					for i in range(nruns, nruns+numprocesses):
						treefile='tree'+str(i)+'.txt'
						alnfile='bootaln'+str(i)+'.fasta'
						p=multiprocessing.Process(target=self.makeAlignmentGT, args=( treefile, prealn_file, alnfile ))
						jobs.append(p)
						p.start()
					for p in jobs:
						p.join()
					nruns+=numprocesses ## now increment
				elif (n-nruns)<numprocesses:
					jobs=[]
					for i in range(n-nruns, n):
						treefile='tree'+str(i)+'.txt'
						alnfile='bootaln'+str(i)+'.fasta'
						p=multiprocessing.Process(target=self.makeAlignmentGT, args=( treefile, prealn_file, alnfile ))
						jobs.append(p)
						p.start()
					for p in jobs:
						p.join()
					break
		pool.terminate()
				
		return 0


class MafftAligner(Aligner):
	def __init__(self, executable, options):
		''' "executable" is the path to the MAFFT and its options are given by "options" '''
		self.executable = executable
		self.options = options
	
	def makeAlignment( self, prealn_file, alnfile):
		print "Making initial alignment with MAFFT"
		align=self.executable+' '+self.options+' '+prealn_file+' > '+alnfile
		runalign=subprocess.call(str(align),shell=True)

		return 0

	def makeAlignmentGT( self, treefile, prealn_file, alnfile):
		align=self.executable+' '+self.options+' --treein '+treefile+' '+prealn_file+' > '+alnfile ## Note that we do not use "--retree 1." Providing an input tree will still capture alignment stochasticity without forcing a poorer alignment, as retree 1 has the potential to do.
		runalign=subprocess.call(str(align), shell=True)
		return 0
			
	def processTrees(self, n, infile):
		''' Takes the bootstrapped trees out from a single file and process each into MAFFT format.'''
		trees=TreeList()
		trees.read(open(infile, 'r'), 'newick', as_rooted=True)
		for i in range(n):
			rawtree = trees[i]
			# Remove any polytomies and update splits since it was forced in as rooted 3 lines above
			rawtree.resolve_polytomies()
			rawtree.update_splits()
			rawtree2=str(rawtree).replace('[&R] ','')
			outtree = "tree"+str(i)+".txt"
			self.Tree2Mafft(rawtree2, outtree)			
		return 0
	
	
	#################################################################################
	### The following two functions are simply the newick2mafft.rb script re-written in python. Note that the MAFFT authors have been notified of this modified code and intend to distribute it.
				
	def killMegatomy(self, tree):
		''' First part of the newick2mafft.rb script.'''
		findMegatomy = re.search(",(\d+):(\d+\.\d*),(\d+):(\d+\.\d*)", tree)
		while findMegatomy:
			hit1 = findMegatomy.group(1)
			hit2 = findMegatomy.group(2)
			hit3 = findMegatomy.group(3)
			hit4 = findMegatomy.group(4)
			wholeMega = ","+hit1+":"+hit2+","+hit3+":"+hit4
			tree = tree.replace(wholeMega, ",XXX")
			poshit = tree.index("XXX")
			i=poshit
			height=0
			while i >= 0:
				if height == 0 and tree[i]=='(':
					break
				if tree[i]==')':
					height+=1
				elif tree[i] == '(':
					height-=1
				i-=1
			poskakko=i
			zenhan=tree[0:poskakko+1]
			treelen = len(tree)
			tree = zenhan + "(" + tree[(poskakko+1):treelen]
			tree=tree.replace('XXX', hit1+":"+hit2+"):0,"+hit3+":"+hit4)
			findMegatomy = re.search(",(\d+):(\d+\.\d*),(\d+):(\d+\.\d*)", tree)	
		return tree


	def Tree2Mafft(self, tree, outfile):	
		''' Second part of the newick2mafft.rb script. Converts a newick tree into MAFFT's native format.'''
		outhandle=open(outfile, 'w')
		
		# Replace forbidden characters. 
		tree = re.sub(" ", "", tree)
		tree = re.sub("\d\.*\d*e-\d+", "0", tree)
		tree = re.sub("\[.*?\]", "", tree)
		tree = re.sub("[_*?:]", ":", tree)
		
		memi=[-1,-1]
		leni=[-1,-1]
		findparen = re.search('\(', tree)
		while findparen:
			tree=self.killMegatomy(tree)
			find = re.search("\((\d+):(\d+\.*\d*),(\d+):(\d+\.*\d*)\)", tree)
			if find:
				memi[0] = find.group(1)
				leni[0] = find.group(2)
				memi[1] = find.group(3)
				leni[1] = find.group(4)
				wholeclade = "("+memi[0]+":"+leni[0]+","+memi[1]+":"+leni[1]+")"
				tree = tree.replace(wholeclade, "XXX")
				if int(memi[1]) < int(memi[0]):
					memi.reverse()
					leni.reverse()
				tree=tree.replace("XXX", memi[0])	
				outhandle.write(memi[0]+' '+memi[1]+' '+leni[0]+' '+leni[1]+'\n')
			findparen = re.search('\(', tree)
		outhandle.close()
		return 0
	##############################################################################################
		
		
		
		

class ClustalAligner(Aligner):
	def __init__(self, executable, options):
		'''"executable" is the path to the CLUSTALW executable, and options are user-specified options.'''
		self.executable = executable
		self.options = options
	
	def makeAlignment(self, prealn_file, alnfile):
		print "Making initial alignment with CLUSTALW"
		align=self.executable+' -align -output=FASTA -outorder=INPUT '+self.options+' -infile='+prealn_file+' -outfile='+alnfile
		runalign=subprocess.call(str(align),shell=True)
	
	def makeAlignmentGT(self, treefile, prealn_file, alnfile):
		align=self.executable+' -align -output=FASTA -outorder=INPUT '+self.options+' -usetree='+treefile+' -infile='+prealn_file+' -outfile='+alnfile
		runalign=subprocess.call(str(align), shell=True)
	
	def processTrees(self, n, infile):
		''' All BS trees need to be taken out of the single file they're in. As different alignment programs have different tree requirements, the formatting can happen in this class'''	
		print "processing trees"
		inhandle=open(infile, 'r')
		lines=inhandle.readlines()
		for i in range(n):
			mytree = str(lines[i])
			outhandle = "tree"+str(i)+".txt"
			outfile = open(outhandle, "w")	
			outfile.write(mytree)	
			outfile.close()	
		return 0



class MuscleAligner(Aligner):
	def __init__(self, executable, options):
		'''"executable" is the path to the MUSCLE executable, and options are user-specified options.'''
		self.executable = executable
		self.options = options
	
	def makeAlignment(self, prealn_file, alnfile):
		print "Making initial alignment with Muscle"
		align=self.executable+' '+self.options+' -in '+prealn_file+' -out temp.aln'
		runalign=subprocess.call(str(align),shell=True)
		self.reorderAlignment('temp.aln', alnfile)
	
	def makeAlignmentGT(self, treefile, prealn_file, alnfile):
		temp = 'temp'+alnfile
		align=self.executable+' '+self.options+' -quiet -usetree_nowarn '+treefile+' -in '+prealn_file+' -out '+temp
		runalign=subprocess.call(str(align), shell=True)
		self.reorderAlignment(temp, alnfile)
	
	def processTrees(self, n, infile):
		''' Places bootstrapped trees into separate files. All trees are rooted here.'''
		trees=TreeList()
		trees.read(open(infile, 'r'), 'newick', as_rooted=True)
		for i in range(n):		
			rawtree = trees[i]
			rawtree.resolve_polytomies()
			rawtree.update_splits()
			mytree=str(rawtree).replace('[&R] ','')
		
			outtree="tree"+str(i)+".txt"
			outfile = open(outtree, 'w')
			outfile.write(mytree+';')
			outfile.close()
			
		return 0
	
	def reorderAlignment(self, infile, outfile):
		'''Given an input alignment in fasta format, reorder the sequences (based on ascending id's which are now all ints. indexing begins at 1.) and rewrite the file, again in fasta. Note that MUSCLE's ability to do this is deprecated, so we have to do it here instead of relying of MUSCLE to do it automatically.'''
		rawaln = AlignIO.read(infile, 'fasta')
		outaln = open(outfile, 'w')
		numseq = len(rawaln)
		for i in range(1, numseq+1):
			for aln in rawaln:
				if aln.id == str(i):
					outaln.write('>'+str(i)+'\n'+str(aln.seq)+'\n')
					break
		outaln.close()
		return 0
				
		
		
		
		
		
		
		
		
		
		
		

		
		
