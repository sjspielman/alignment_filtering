import re, os, sys, subprocess, fnmatch
import shutil
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment


def GetFiles(ext, dirpath):
	files=[]
	dir=os.listdir(dirpath)
	for file in dir:
		if fnmatch.fnmatch(file, '*.'+str(ext)):
			files.append(file)
	return files
########

def mytranslate(seq, codon_table):
	peptide = ''
	for i in xrange(0, len(seq), 3):
		amino_acid=''
		codon = seq[i: i+3]
		if codon=='---':
			amino_acid='-'
		elif codon=='NNN':
			amino_acid='X'
		else:
		 	amino_acid = codon_table.get(codon, '')
		if amino_acid != '':
			peptide += amino_acid
		else:
			print codon
			print "FAAAAIIIILLLL"
			break
	return peptide
    
bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))
## codon_table['ATG'], you'll get M 

rundict={1:'med128b', 2:'med16p', 3:'med64r', 4:'short16b', 5:'short64p', 6:'tiny16b', 7:'tiny64p', 8:'med128p', 9:'med16r', 10:'short128b', 11:'short16p', 12:'short64r', 13:'tiny16p', 14:'tiny64r', 15:'med128r', 16:'med64b', 17:'short128p', 18:'short16r', 19:'tiny128p', 20:'tiny16r', 21:'med16b', 22:'med64p', 23:'short128r', 24:'short64b', 25:'tiny128r', 26:'tiny64b', 27:'tiny128b'}

for n in range(1,28):

	basedir='guided_'+rundict[n]

	newdir='aa_guided_'+rundict[n]
	os.mkdir(newdir)

	files=GetFiles('phy',basedir)
	for file in files:
		print file
		find=re.search('([_\w\d]+\d+)\.phy', file)
		if find:
			outfile=newdir+'/'+find.group(1)+'.phy'
			#Translate each file and save to new dir
			inhandle=open(basedir+'/'+file, 'r')
			parsed=AlignIO.read(inhandle, 'phylip')
			inhandle.close()
			newMSA=MultipleSeqAlignment([])
			for record in parsed:
				newseq = mytranslate(str(record.seq), codon_table)
				record=SeqRecord(Seq(str(newseq), generic_protein), id=record.id, description='')
				newMSA.append(record)
			outhandle=open(outfile, 'w')
			umm=AlignIO.write(newMSA, outhandle, 'phylip')
			outhandle.close()
	
	