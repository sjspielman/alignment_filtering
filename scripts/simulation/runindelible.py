# 8/7/13. Simulate sequences with Indelible.

import re, os, sys, subprocess, shutil
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,generic_protein
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment


def makeDir(dir):
	if not os.path.exists(dir):
		os.mkdir(dir)
	
##### Translates and allows for gaps
def myTranslate(seq, codon_table):
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
### Relevant to the myTranslate function.
bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))



simdir = '/Users/sjspielman/Dropbox/aln/results/Simulation/'
rundir = os.getcwd()

type = sys.argv[1] # neutral or HA
if type == 'neutral':
	outdir_base = simdir+'Neutral/'
elif type == 'HA':
	outdir_base = simdir+'HA/'
else:
	assert(1==0), "Need to specify argument neutral or HA"
makeDir(outdir_base)
	
	

# Simulate according to HA params, but with div/indels of the trees and their original alignments.
#gene_indels={'or5':'0.053', 'rho':'0.019', 'prk':'0.0041', 'flat':'0.0066'}
gene_indels = {'rho':0.019}

for gene in gene_indels:
	cfile=open(type+'_control.txt', 'r')
	control=cfile.read()
	cfile.close()
		
	makeDir(outdir_base+'sequences/')
	makeDir(outdir_base+'truerates/')
	makeDir(outdir_base+'logs/')
	
	seqbase=outdir_base+'sequences/'+gene+'/'
	ratesbase=outdir_base+'truerates/'+gene+'/'
	logbase=outdir_base+'logs/'+gene+'/'
	
	makeDir(seqbase)
	makeDir(ratesbase)
	makeDir(logbase)
	
	seqdir=seqbase+'seqs/'
	ratesdir=ratesbase+'seqs/'
	logdir = logbase+'seqs/'
	
	makeDir(seqdir)
	makeDir(ratesdir)
	makeDir(logdir)
	
	treefile = simdir+'GenerateTrees/real_'+gene+'.tre'
	infile = open(treefile, 'r')
	tree = infile.read()
	infile.close()
			
	## Find and remove any scientific notation from the tree
	tree=tree.rstrip()
	tree = re.sub('\d\.\d+e-\d+', '0', tree) # just replace it with a 0. if it had to resort to scientific notation, then it's fine.
	
	control=re.sub(r'\[TREE\].+;', r'[TREE] mytree '+tree, control)
				
	cfile=open('control.txt', 'w')
	cfile.write(control)
	cfile.close()
		
	for n in range(100): ## 100 reps for each
		os.chdir(rundir)
		runIndelible = subprocess.call('indelible control.txt', shell=True)

		# New file names
		rawsim_codon='rawsim_codon'+str(n)+'.fasta'
		rawsim_aa='rawsim_aa'+str(n)+'.fasta'
		outtrue_nuc='truealn_codon'+str(n)+'.fasta'
		outtrue_aa='truealn_aa'+str(n)+'.fasta'
		ratefile='truerates'+str(n)+'.txt'
	
		# Copy output to correct directory and rename
		shutil.move('results_TRUE.phy', seqdir)
		shutil.move('results.fas', seqdir+rawsim_codon)
		shutil.move('results_RATES.txt', ratesdir+ratefile)
		shutil.move('LOG.txt', logdir+'log'+str(n)+'.txt')
		
		os.chdir(seqdir)
	
		#Translate rawsim file to be "rawsim_aa.fasta" 
		rawtrans=open(rawsim_aa, 'w')
		rawhandle=open(rawsim_codon, 'r')
		rawparsed=list(SeqIO.parse(rawhandle, 'fasta'))
		rawhandle.close()
		for record in rawparsed:
			numgap = str(record.seq).count('-')
			assert(numgap != len(record.seq)), "Too many gaps. Indel rate is poorly chosen."
			newseq = record.seq.translate()
			rawtrans.write('>'+record.id+'\n'+str(newseq)+'\n')
		rawtrans.close()
		
		
		# Fix the truealn formats. Save both a truealn_codon.fasta and a truealn_aa.fasta
		inhandle=open('results_TRUE.phy', 'r')
		lines=inhandle.readlines()
		inhandle.close()
		newMSA_nuc=MultipleSeqAlignment([])
		newMSA_aa=MultipleSeqAlignment([])	
		
		for i in range(1, len(lines)):   ##skip the first line since it's all phylip-y
			find=re.search('^(\w+)\s+([\w-]+)', lines[i])
			if find:
				id=find.group(1)
				seq_nuc=find.group(2)
				seq_aa=myTranslate(seq_nuc, codon_table)
				record_nuc=SeqRecord(Seq(str(seq_nuc), generic_dna), id=str(i), description='')
				newMSA_nuc.append(record_nuc)
				record_aa=SeqRecord(Seq(str(seq_aa), generic_dna), id=str(i), description='')
				newMSA_aa.append(record_aa)
		outhandle=open(outtrue_nuc, 'w')
		umm=AlignIO.write(newMSA_nuc, outhandle, 'fasta')
		outhandle.close()
		outhandle=open(outtrue_aa, 'w')
		umm=AlignIO.write(newMSA_aa, outhandle, 'fasta')
		outhandle.close()
		
		os.remove('results_TRUE.phy')

os.chdir(rundir)	
os.remove('control.txt')
os.remove('trees.txt')
