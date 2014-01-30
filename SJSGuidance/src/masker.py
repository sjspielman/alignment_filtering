from numpy import *
from Bio import AlignIO, SeqIO
from Bio.Alphabet import *
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
import re

class Masker:
	def __init__(self, bootstrapper):
		'''initialization function'''
		self.bootstrapper=bootstrapper
		return			
	
	def notGapSites(self, parsed, col_index):
		## Determine which sites in a column are gaps. Returns a list of lists whereby each nested list contains the indices for which ones are not gaps or a 'G' for if that position is a gap
		notGaps=[]
		findgaps = str(parsed[:,col_index])
		for counter in range(len(findgaps)): #column!
			if findgaps[counter]!='-':
				notGaps.append(counter)
			else:
				notGaps.append('G')
		return notGaps


	def maskResidues(self, refMSA_file, numseq, alnlen, scores, x, formatout, final_file, seqType):
		''' Masks poorly aligned residues whose score is <x. Will NOT mask gaps.'''
		
		new='?'	
		parsed = AlignIO.read(refMSA_file, 'fasta')
		while fracmasked > maxmasked:
			newseqs=[]
			numres=0
			totalmasked=0
			maskedMSA=MultipleSeqAlignment([])
			for row in range(numseq):
				newseq=''
				for position in range(alnlen):
					thispos=str(parsed[row].seq[position])
					isgap=re.search('-', thispos)
					if isgap:
						newseq=newseq+parsed[row].seq[position]
					else:
						numres+=1
						thescore=scores[row][position]
						if thescore<x: #mask if shitty
							newseq=newseq+new
							totalmasked+=1
						else: #or, keep that position
							newseq=newseq+parsed[row].seq[position]
				newseqs.append(newseq)
		
		for i in range(numseq):
			if str(seqType)=='protein':
				aln_record=SeqRecord(Seq(newseqs[i],generic_protein), id=str(i+1), description='')
			elif str(seqType)=='nucleotide':
				aln_record=SeqRecord(Seq(newseqs[i],generic_dna),id=str(i+1), description='')
			maskedMSA.append(aln_record)
	
		outhandle=open(final_file, 'w')
		outhandle.write(maskedMSA.format(str(formatout)))
		outhandle.close()
		
		return (totalmasked)	
