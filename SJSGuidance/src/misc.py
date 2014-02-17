#!/usr/bin/python

# This file contains miscellaneous functions used throughout SJSGuidance.

from Bio import SeqIO
from Bio import AlignIO, SeqIO
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

from numpy import *
import os
import subprocess

def prepareDir(directory, save=False, newname=None):
	''' Double check that BootDir exists and clear it of all files, as needed. ''' 
	if save:
		# Be sure newname is properly named. Yes, not an ideal way. I know. But it's fine.
		if ".tgz" not in newname:
			newname+=".tgz"
		compress = "tar -czf "+newname+" "+directory
		subprocess.call(compress, shell=True)
	
	# Cleanup and/or make
	if (os.path.exists(directory)):
		bootfiles=os.listdir(directory)
		for file in bootfiles:
			os.remove(directory+file)	
	else:
		os.mkdir(directory)
		
		
		
		
		
def buildMap(unaligned_aa, format, prealn_file):	
	''' Create map to original ids. New ids are ints'''	
	idmap={}		
	count=0
	
	infile=open(unaligned_aa, 'rU')
	parsed = list(SeqIO.parse(infile, str(format)))
	infile.close()		
	
	out=open(prealn_file, 'w')
	for record in parsed:	
		count+=1
		idmap[count] = str(record.id)
		seq=str(record.seq)
		out.write('>'+str(count)+'\n'+seq+'\n')
	out.close()
	return idmap

def unMap(idmap, alnfile, final_alnfile, numseq):
	''' Write alignment file with correct ids '''
	
	aln = AlignIO.read(alnfile, "fasta")
	out = open(final_alnfile, 'w')
	for i in range(numseq):
		out.write(">"+idmap[i+1]+"\n"+str(aln[i].seq)+"\n")
	out.close()


def maskResidues(refMSA_file, numseq, alnlen, map, scores, x, formatout, final_file, seqType):
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
			if thispos=='-':
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
		if str(seqType)=='prot':
			aln_record=SeqRecord(Seq(newseqs[i],generic_protein), id=str(map[i+1]), description='')
		elif str(seqType)=='dna':
			aln_record=SeqRecord(Seq(newseqs[i],generic_dna), id=str(map[i+1]), description='')
		maskedMSA.append(aln_record)

	outhandle=open(final_file, 'w')
	outhandle.write(maskedMSA.format(str(formatout)))
	outhandle.close()





def Pal2Nal(palfile, nucfile, paltype, nuctype, outfile, outputformat):
	''' Convert a protein alignment to a nucleotide alignment. Can handle the ambiguities N, X, ? '''
	## Arguments:
	##	-> palfile = protein alignment file
	##  -> nucfile = unaligned nucleotide file. Sequences should have same name and be in same order as in palfile
	##  -> paltype = palfile format (eg fasta...)
	##  -> nuctype = nucfile format
	##  -> outfile = output file for nucleotide alignment
	##  -> outputformat = format for output nucleotide alignment
	
	
	aln_parsed=list(SeqIO.parse(str(palfile), str(paltype)))
	nuc_parsed=list(SeqIO.parse(str(nucfile), str(nuctype)))	
	
	
	assert ( len(aln_parsed) == len(nuc_parsed) ), "The protein and nucleotide files have a different number of sequences. Are you sure they correspond?"
	numseq=len(aln_parsed)
	
	nucMSA=MultipleSeqAlignment([])
	for p in range(0,numseq):
		pal_seq=str(aln_parsed[p].seq) #aa alignment sequence
		pal_id=str(aln_parsed[p].id)
		for n in range(0,numseq):
			if nuc_parsed[n].id==pal_id:
				nuc_seq=str(nuc_parsed[n].seq)
				nal=str()
				start=0 #counter for codon starting position
				end=3 #counter for codon ending position
				for position in pal_seq:
					#If gapped, missing, or masked position in alignment, append 3 gaps/missing/NNN to new string 
					if position=='-':
						codon='---'
						nal=nal+codon
					elif position=='?':
						codon='???'
						nal=nal+codon
						start+=3
						end+=3
					elif position=='X':
						codon='NNN'
						nal=nal+codon
						start+=3
						end+=3
					#If amino acid there, append corresponding codon
					else:
						codon=str(nuc_seq[start:end])
						nal=nal+codon
						start+=3
						end+=3
				#Make nucleotide MSA object
				nal_seq=Seq(nal)
				aln_record=SeqRecord(Seq(nal, generic_dna), id=pal_id, description='')
				nucMSA.append(aln_record)
			else:
				continue
				
	#write alignment to file
	outfile=open(outfile, 'w')
	umm=AlignIO.write(nucMSA, outfile, "fasta")
	outfile.close()
	return 0



	