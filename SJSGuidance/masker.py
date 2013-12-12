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


	def maskResidues(self, refMSA_file, numseq, alnlen, scores, x, idmap, formatout, final_file, seqType, blah, save_x_file, alg):
		''' Masks poorly aligned residues whose score is <x. Will NOT mask gaps.'''
		
		maxmasked = 0.4
		
		## We want to be sure that the roundings are getting there, so if x=0.9, we want all the >=0.895. Thus, mask at x-0.005 !!
		x = float(x - 0.005)	
		
		fracmasked = 1
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
			
			fracmasked = float(totalmasked)/float(numres)
			if fracmasked <= maxmasked:
				break
			if fracmasked <= (maxmasked+0.05) and fracmasked > maxmasked:
				x-=0.05
			elif fracmasked > (maxmasked+0.05):
				x-=0.1

		
		for i in range(numseq):
			if str(seqType)=='protein':
				aln_record=SeqRecord(Seq(newseqs[i],generic_protein), id=str(idmap[i+1]), description='')
			elif str(seqType)=='nucleotide':
				aln_record=SeqRecord(Seq(newseqs[i],generic_dna),id=str(idmap[i+1]), description='')
			maskedMSA.append(aln_record)
	
		outhandle=open(final_file, 'w')
		outhandle.write(maskedMSA.format(str(formatout)))
		outhandle.close()
		
		xfile=open(save_x_file, 'a')
		xfile.write(str(blah)+'\t'+str(alg)+'\t'+str(x)+'\t'+str(fracmasked)+'\n')
		xfile.close()
		
		return (totalmasked)	

	
	def maskColumns(self, refMSA_file, numseq, alnlen, scores, x, idmap, formatout, final_file, seqType):
		''' Masks poorly aligned columns whose score <x. Unweighted average over residue scores (weighting doesn't change anything more than 1%).'''
		
		#All columns masked to become ?
		mcol=numseq*'?'
		x=float(x)
		
		parsed = AlignIO.read(refMSA_file, 'fasta')
		raw=[]
		
		## Loop over columns in the scores.		
		col_counter=0
		for pos_score in scores.T: #transposes it so loops over columns of all_scores instead of default rows
		
			notGaps=self.notGapSites(parsed, col_counter) # Making a bunch of lists containing indices for if it's not a gap, 'G' if it is a gap. 
			count=0
			sum=0
			for i in range(len(notGaps)):
				if notGaps[i] != 'G':
					sum+=pos_score[i]			
					count+=1
			
			colscore=float(sum)/float(count)				
			
			colseq=''	
			for row in range(numseq):
				colseq=colseq+parsed[row].seq[col_counter]
								
			if colscore>=x:
				raw.append(str(colseq))
			elif colscore<x:
				raw.append(str(mcol))			
			col_counter+=1
		
		## Here, I turn them into lists of list. Each list is a ROW in the final MSA. hurray!
		raw_transposed=zip(*raw)


		## The alignments have been masked. Now need to unmap the names.
		maskedMSA=MultipleSeqAlignment([])
		map_counter=0
		for row in raw_transposed:
			newrow=''.join(row)
			if str(seqType)=='protein':
				aln_record=SeqRecord(Seq(newrow,generic_protein),id=str(idmap[map_counter+1]), description='')
			elif str(seqType)=='nucleotide':
				aln_record=SeqRecord(Seq(newrow,generic_dna),id=str(idmap[map_counter+1]), description='')
			maskedMSA.append(aln_record)
			map_counter+=1

		
		outhandle=open(final_file, 'w')
		outhandle.write(maskedMSA.format(str(formatout)))
		outhandle.close()
		return 0
