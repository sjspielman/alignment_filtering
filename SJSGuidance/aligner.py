import subprocess, sys, re, os, multiprocessing
from multiprocessing import pool
from dendropy import *
from Bio import AlignIO
    
class Aligner:
	def __init__(self):
		'''initialization function'''
		return

	def makeAlignment( self, prealn_file, alnfile ):
		'''makes a standard alignment'''
		print "Aligner base class. Function makeAlignment() is not implemented."

	def makeAlignmentGT( self, treefile, prealn_file, alnfile):
		'''makes an alignment using the provided guide tree'''
		print "Aligner base class. Function makeAlignmentGT() is not implemented."
	
	
	
	def multiMakeAlignmentsGT(self, prealn_file, n, numprocesses):
		'''This function creates a pool of processes and then executes them on one fewer CPU cores than the machine has.'''
		if numprocesses=='':
			numprocesses=multiprocessing.cpu_count() - 1 # leave a spare, you never know..		
		jobs=[]		
		pool=multiprocessing.Pool(numprocesses) ## can also for the argument here do multiprocessing.cpu_count()
		
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
		'''"executable" is the path to the MAFFT executable and options are given by options'''
		self.executable = executable
		self.options = options
	
	def makeAlignment( self, prealn_file, alnfile):
		'''makes a standard alignment'''
		print "Making initial alignment with MAFFT"
		align=self.executable+' '+self.options+' '+prealn_file+' > '+alnfile
		runalign=subprocess.call(str(align),shell=True)

		return 0

	def makeAlignmentGT( self, treefile, prealn_file, alnfile):
		'''makes an alignment using the provided guide tree'''
		#print "Making alignment with MAFFT from guide tree, using ", self.executable
		align=self.executable+' '+self.options+' --retree 1 --treein '+treefile+' '+prealn_file+' > '+alnfile ## TAUguidance forces retree 1, which is some more academic dishonesty, but what else is new.
		runalign=subprocess.call(str(align), shell=True)
		return 0
			
	def processTrees(self, n, infile):
		''' All BS trees need to be taken out of the single file they're in. As different alignment programs have different tree requirements, the formatting can happen in this class'''
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
	
	### THIS FUNCTION WILL CALL THE RUBY SCRIPT. IT IS NOT USED AS THE RUBY SCRIPT IS CODED HERE IN PYTHON (line for line). KEEP IN CASE THOUGH!! ###
	def processTreesRUBY(self, n, infile):
		'''NOT USED but just in case, it's here.'''
		trees=TreeList()
		trees.read(open(infile, 'r'), 'newick', as_rooted=True)
		for i in range(n):
			rawtree = trees[i]
			# Remove any polytomies and assure that it is rooted.
			###rawtree.resolve_polytomies() # probably not necessary, but just in case.
			###rawtree.reroot_at_midpoint()
			rawtree.update_splits()
			rawtree2=str(rawtree).replace('[&R] ','')
			
			temptree = open('temp.txt', 'w')
			temptree.write(rawtree2+';')
			temptree.close()
			
			outtree = "tree"+str(i)+".txt"
			ruby = subprocess.call('ruby ../newick2mafft.rb temp.txt > '+outtree, shell=True)			
		return 0
		
	
		
	def killMegatomy(self, tree):
		''' first chunk of that ruby script'''
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
		''' converts a newick tree into mafft's native format. line-by-line the ruby script. We have sent to Katoh too, so now mafft will release this script.'''
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
		
		
		
		
		

class ClustalAligner(Aligner):
	def __init__(self, executable, options):
		'''"executable" is the path to the CLUSTALW executable'''
		self.executable = executable
		self.options = options
	
	def makeAlignment(self, prealn_file, alnfile):
		'''makes a standard alignment'''
		print "Making initial alignment with CLUSTALW"
		align=self.executable+' -align -output=FASTA -outorder=INPUT '+self.options+' -infile='+prealn_file+' -outfile='+alnfile
		runalign=subprocess.call(str(align),shell=True)
	
	def makeAlignmentGT(self, treefile, prealn_file, alnfile):
		'''makes an alignment using the provided guide tree'''
		align=self.executable+' -align -output=FASTA -outorder=INPUT '+self.options+' -usetree='+treefile+' -infile='+prealn_file+' -outfile='+alnfile
		runalign=subprocess.call(str(align), shell=True)
	
	def processTrees(self, n, infile):
		''' All BS trees need to be taken out of the single file they're in. As different alignment programs have different tree requirements, the formatting can happen in this class'''	
		print "processing trees"
		inhandle=open(infile, 'r')
		lines=inhandle.readlines()
		for i in range(n):
			mytree = str(lines[i])
			print mytree
			outhandle = "tree"+str(i)+".txt"
			outfile = open(outhandle, "w")	
			outfile.write(mytree)	
			outfile.close()	
		return 0



class MuscleAligner(Aligner):
	def __init__(self, executable, options):
		'''"executable" is the path to the CLUSTALW executable'''
		self.executable = executable
		self.options = options
	
	def makeAlignment(self, prealn_file, alnfile):
		'''makes a standard alignment'''
		print "Making initial alignment with Muscle"
		align=self.executable+' '+self.options+' -in '+prealn_file+' -out temp.aln'
		runalign=subprocess.call(str(align),shell=True)
		self.reorderAlignment('temp.aln', alnfile)
	
	def makeAlignmentGT(self, treefile, prealn_file, alnfile):
		'''makes an alignment using the provided guide tree'''
		temp = 'temp'+alnfile
		align=self.executable+' '+self.options+' -quiet -usetree_nowarn '+treefile+' -in '+prealn_file+' -out '+temp
		runalign=subprocess.call(str(align), shell=True)
		self.reorderAlignment(temp, alnfile)
	
	def processTrees(self, n, infile):
		''' All BS trees need to be taken out of the single file they're in and rooted.'''
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
		'''Given an input alignment in fasta format, reorder the sequences (based on ascending id's which are now all ints. indexing begins at 1.) and rewrite the file, again in fasta.'''
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
				

