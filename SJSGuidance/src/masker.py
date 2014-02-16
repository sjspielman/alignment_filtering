from numpy import *
from Bio import AlignIO, SeqIO
from Bio.Alphabet import *
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
import re

def maskResidues(self, refMSA_file, numseq, alnlen, scores, x, formatout, final_file, seqType):
	''' Masks poorly aligned residues whose score is <x. Will NOT mask gaps.'''
	
	new='?'	
	parsed = AlignIO.read(refMSA_file, 'fasta')
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
				if round(thescore)<x: #mask if below threshold. use round to ensure we get the ones >= 0.5 below threshold
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
