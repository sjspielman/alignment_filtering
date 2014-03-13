import subprocess, sys, re, os, multiprocessing
from dendropy import *
from Bio import AlignIO
    
class Aligner:
	def __init__(self, **kwargs):
		self.prealn_file = kwargs.get("prealn_file", "prealn.fasta")
		self.refaln_file = kwargs.get("refaln_file", "refaln.fasta")
		self.n           = kwargs.get("bootstraps", 100)
		self.cpu         = kwargs.get("cpu", 1) 

	
	def multiMakeAlignmentsGT(self):
		'''Makes n bootstrap alignments. '''
		# Note that this function is rather clunky because multiprocessing.pool() does not work when inside a class in python. The code here is a decent(?) workaround of this unfortunate issue.


		pool=multiprocessing.Pool(self.cpu)
		jobs=[]
		
		if self.cpu >= self.n:
			for i in range( self.n ):
				treefile='tree'+str(i)+'.txt'
				alnfile='bootaln'+str(i)+'.fasta'
				p=multiprocessing.Process(target=self.makeAlignmentGT, args=( treefile, self.prealn_file, alnfile ))
				jobs.append(p)
				p.start()
			for p in jobs:
				p.join()
		
		elif self.cpu < self.n:
			nruns=0
			while ( nruns < self.n ):
				if ( self.n - nruns ) >= self.cpu:
					jobs=[]
					for i in range( nruns, nruns + self.cpu ):
						treefile='tree'+str(i)+'.txt'
						alnfile='bootaln'+str(i)+'.fasta'
						p=multiprocessing.Process(target=self.makeAlignmentGT, args=( treefile, self.prealn_file, alnfile ))
						jobs.append(p)
						p.start()
					for p in jobs:
						p.join()
					nruns += self.cpu
				elif ( self.n - nruns ) < self.cpu:
					jobs=[]
					for i in range( self.n - nruns, self.n ):
						treefile='tree'+str(i)+'.txt'
						alnfile='bootaln'+str(i)+'.fasta'
						p=multiprocessing.Process(target=self.makeAlignmentGT, args=( treefile, self.prealn_file, alnfile ))
						jobs.append(p)
						p.start()
					for p in jobs:
						p.join()
					break
		pool.terminate()



	
class MafftAligner(Aligner):
	''' Implement alignments with mafft '''
	def __init__(self, **kwargs):
		super(MafftAligner, self).__init__(**kwargs)
		self.executable = kwargs.get("executable", "mafft")
		self.options    = kwargs.get("options", " --auto --quiet ")
	
	def makeAlignment( self, prealn_file, alnfile):
		print "Making initial alignment with MAFFT"
		align=self.executable+' '+self.options+' '+self.prealn_file+' > '+self.refaln_file
		runalign=subprocess.call(str(align),shell=True)


	def makeAlignmentGT( self, treefile, prealn_file, alnfile):
		align=self.executable+' '+self.options+' --treein '+treefile+' '+self.prealn_file+' > '+alnfile ## Note that we do not use "--retree 1." Providing an input tree will still capture alignment stochasticity without forcing a poorer alignment, as retree 1 has the potential to do.
		runalign=subprocess.call(str(align), shell=True)

			
	def processTrees(self, infile):
		''' Takes the bootstrapped trees out from a single file and process each into MAFFT format.'''
		trees=TreeList()
		trees.read(open(infile, 'r'), 'newick', as_rooted=True)
		for i in range( self.n ):
			rawtree = trees[i]
			# Remove any polytomies and update splits since it was forced in as rooted 3 lines above
			rawtree.resolve_polytomies()
			rawtree.update_splits()
			rawtree2=str(rawtree).replace('[&R] ','')
			outtree = "tree"+str(i)+".txt"
			self.Tree2Mafft(rawtree2, outtree)			

	
	
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
		
		
		

########################### You are welcome to run this software with either of these two classes, although we don't recommend. Alignment made with "mafft --auto" are much better than muscle or clustal alignments!! ###############
class ClustalAligner(Aligner):
	''' Implement alignments with clustalw '''
	def __init__(self, **kwargs):
		super(ClustalAligner, self).__init__(**kwargs)
		self.executable = kwargs.get("executable", "clustalw2")
		self.user_options    = kwargs.get("options", " -quiet ")
		self.mand_options = ' -align -output=FASTA -outorder=INPUT ' 

	
	def makeAlignment(self, prealn_file, alnfile):
		print "Making initial alignment with CLUSTALW"
		align=self.executable + self.mand_options + self.user_options + ' -infile=' + self.prealn_file + ' -outfile=' + self.refaln_file
		runalign=subprocess.call(str(align),shell=True)
	
	def makeAlignmentGT(self, treefile, alnfile):
		align=self.executable + self.mand_options + self.user_options + ' -usetree=' + treefile + ' -infile=' + self.prealn_file + ' -outfile=' + alnfile
		runalign=subprocess.call(str(align), shell=True)
	
	def processTrees(self, infile):
		''' All bootstrap trees need to be taken out of the single file they're in. As different alignment programs have different tree requirements, the formatting can happen in this class'''	
		print "Processing bootstrap trees"
		inhandle=open(infile, 'r')
		lines=inhandle.readlines()
		for i in range( self.n ):
			mytree = str(lines[i])
			outhandle = "tree"+str(i)+".txt"
			outfile = open(outhandle, "w")	
			outfile.write(mytree)	
			outfile.close()	
		return 0



class MuscleAligner(Aligner):
	def __init__(self, **kwargs):
		super(MuscleAligner, self).__init__(**kwargs)
		self.executable = kwargs.get("executable", "muscle")
		self.user_options    = kwargs.get("options", " -quiet ")
	
	def makeAlignment(self, prealn_file, alnfile):
		print "Making initial alignment with Muscle"
		align=self.executable+' '+self.options+' -in '+prealn_file+' -out temp.aln'
		runalign=subprocess.call(str(align),shell=True)
		self.reorderAlignment('temp.aln', alnfile)
	
	def makeAlignmentGT(self, treefile, alnfile):
		temp = 'temp'+alnfile
		align=self.executable+' '+self.options+' -usetree_nowarn '+treefile+' -in '+self.prealn_file+' -out '+temp
		runalign=subprocess.call(str(align), shell=True)
		self.reorderAlignment(temp, alnfile)
	
	def processTrees(self, infile):
		''' Places bootstrapped trees into separate files. All trees are rooted here.'''
		trees=TreeList()
		trees.read(open(infile, 'r'), 'newick', as_rooted=True)
		for i in range(self.n):		
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
				
		
		
		
		
		
		
		
		
		
		
		

		
		
