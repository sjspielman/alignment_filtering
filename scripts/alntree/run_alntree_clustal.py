### 6/12/13
# Masking 70,90 Only the guidance+gweights algorithms
# Uses clustalw for all things.

import re, os, sys, subprocess, shutil
from aligner import *
from treebuilder import *
from scorer import *
from bootstrapper import *
from map import *
from treebuilder import *
from numpy import *
from dendropy import *
from masker import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,generic_protein
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment


################################################################################################
################################################################################################
def Pal2Nal(palfile, nucfile, paltype, nuctype, outfile, outputformat):

	aln_parsed=list(SeqIO.parse(str(palfile), str(paltype)))
	nuc_parsed=list(SeqIO.parse(open(nucfile, 'rU'), str(nuctype)))
	
	if len(aln_parsed)!=len(nuc_parsed):
		print palfile+' '+nucfile+' have different number of sequences! Please make sure that these two files correspond and that all stop codons are removed.'
		assert 1==0
	else:
		numseq=len(aln_parsed)
	
	nucMSA=MultipleSeqAlignment([])
	for p in range(0,numseq):
		pal_seq=str(aln_parsed[p].seq) #aa alignment sequence
		pal_id=str(aln_parsed[p].id)
		for n in range(0,numseq):
			if nuc_parsed[n].id==pal_id:
				nuc_seq=str(nuc_parsed[n].seq)
		print len(pal_seq)
		print pal_seq
		nal=str()
		start=0 #counter for codon starting position
		end=3 #counter for codon ending position
		for position in pal_seq:
			#If gapped, missing, or masked position in alignment, append 3 gaps/missing/NNN to new string 
			if position=='-':
				nal=nal+'---'
			elif position=='?':
				nal=nal+'???'
			elif position=='X':
				nal=nal+'NNN'
			#If amino acid there, append corresponding codon
			else:
				codon=str(nuc_seq[start:end])
				nal=nal+codon
				start+=3
				end+=3
		#Make nucleotide MSA object
		aln_record=SeqRecord(Seq(str(nal), generic_dna), id=pal_id, description='')
		print len(str(nal))
		print aln_record
		print str(nal)
		print '\n\n'
		nucMSA.append(aln_record)
	#write alignment to file
	outfile=open(outfile, 'w')
	umm=AlignIO.write(nucMSA, outfile, str(outputformat))
	outfile.close()
	return 0
	

################################################################################################
################################################################################################


## Keep first blank for indexing with array jobs.
dir=['', 'bigindel_16rc', 'medindel_16rc', 'shortindel_16rc', 'bigindel_64rc', 'medindel_64rc', 'shortindel_64rc', 'bigindel_16p', 'medindel_16p', 'shortindel_16p', 'bigindel_64p', 'medindel_64p', 'shortindel_64p', 'bigindel_16r', 'medindel_16r', 'shortindel_16r', 'bigindel_64r', 'medindel_64r', 'shortindel_64r', 'bigindel_16b', 'medindel_16b', 'shortindel_16b', 'bigindel_64b', 'medindel_64b', 'shortindel_64b']

direc=dir[int(sys.argv[1])]
command='cp -r /home/sjs3495/current/rawdata_5.7/'+str(direc)+' .'
call=subprocess.call(command, shell=True)
alndir_aa='aaguided_'+direc
alndir_nuc='nucguided_'+direc
treedir='aatrees_'+direc
os.mkdir(alndir_aa)
os.mkdir(alndir_nuc)
os.mkdir(treedir)



for blah in range(100):

	os.chdir(direc)
	
	raw='rawsim_aa'+str(blah)+'.fasta'
	rawnuc='rawsim_codon'+str(blah)+'.fasta'
	trueaa='truealn_aa'+str(blah)+'.fasta'
	truenuc='truealn_codon'+str(blah)+'.fasta'
	shutil.copy(trueaa, '../'+alndir_aa)
	shutil.copy(truenuc, '../'+alndir_nuc)
	shutil.copy(raw, '../')
	shutil.copy(raw, '../BootDir')
	shutil.copy(rawnuc, '../BootDir') #need to also bring this into BootDir since will Pal2Nal at the very end there.
	os.chdir('../')				


	#Guidance files
	numprocesses=10
	n = 100  # bootstrap n times
	prealn_file='prealn.fasta'
	refaln_file='refaln.fasta'
	weightfile='treeweights.txt'
	scoreTree_file='scoringtree.tre'
	BootDir='BootDir/'
	#Guidance modules
	amod = CLUSTALAligner("/home/sjs3495/bin/clustalw2/clustalw2", "-quiet -output=FASTA" )
	tmod = builderFastTree( "/share/apps/fasttree-2.1.3/FastTree", " -n "+str(n)+" -noml -nopr -quiet -fastest -nosupport", "-fastest -quiet -nosupport") #last argument = options for scoring tree
	
	mapmod=Map()
	smod = Scorer()
	bmod = Bootstrapper(amod, tmod, smod)
	mmod = Masker(bmod)

	#Final output file names, mask residues only. gaps are maintained. as shown by short64r on 4/28/13, no difference.
	
	#Final output files, 90% masking
	finalr_guidance90="res90_guidance"+str(blah)+".fasta"
	finalr_gweights90="res90_gweights"+str(blah)+".fasta"
	finalc_guidance90="col90_guidance"+str(blah)+".fasta"	
	finalc_gweights90="col90_gweights"+str(blah)+".fasta"	
	
	#Final output files, 70% masking
	finalr_guidance70="res70_guidance"+str(blah)+".fasta"
	finalr_gweights70="res70_gweights"+str(blah)+".fasta"
	finalc_guidance70="col70_guidance"+str(blah)+".fasta"	
	finalc_gweights70="col70_gweights"+str(blah)+".fasta"	
	
	
	
	##############################################################################################################
	##################################### RUNNING "GUIDANCE" HERE ################################################
	print "Starting the Guidances"
	
	#Create initial map for taxa names -> ints
	map=mapmod.ids2int(raw, 'fasta', prealn_file)	

	#Build initial MSA
	amod.makeAlignment(prealn_file, refaln_file)
	
	#Build scoring tree and get weights. Returns dictionary where weights are in taxon order.
	ordered_weights = tmod.buildScoreTree(refaln_file, scoreTree_file, weightfile)

	# Go into BootDir and take necessary files with me
	shutil.copy(refaln_file, BootDir)
	shutil.copy(prealn_file, BootDir)
	shutil.copy(weightfile, BootDir)
	os.chdir(BootDir)

	#SCORING.
	(numseq, alnlen, gscores, gweightscores)=bmod.bootstrap(prealn_file, refaln_file, n, weightfile, ordered_weights, numprocesses)	
	
	
	
	###### MASKING SCORES HERE ########################
	print "Masking residues"
	
	
	temp_res='tempaln_res.aln'
	temp_col='tempaln_col.aln'
	
	
	
	#######Guidance original
	
	## Mask at 90
	
	#residue masking
	mmod.maskResidues(refaln_file, numseq, alnlen, gscores, 0.9, map, 'fasta', temp_res, "protein")
	Pal2Nal(temp_res, rawnuc, 'fasta', 'fasta', finalr_guidance90, 'fasta')
	shutil.copy(finalr_guidance90, '../'+alndir_nuc)
	shutil.copy(temp_res, '../'+alndir_aa+'/'+finalr_guidance90)
	
	#column masking. 
	mmod.maskColumns(refaln_file, numseq, alnlen, gscores, 0.9, map, 'fasta', temp_col, "protein")
	Pal2Nal(temp_col, rawnuc, 'fasta', 'fasta', finalc_guidance90, 'fasta')
	shutil.copy(finalc_guidance90, '../'+alndir_nuc)
	shutil.copy(temp_col, '../'+alndir_aa+'/'+finalc_guidance90)
	
	##Mask at 70
	
	#residue masking
	mmod.maskResidues(refaln_file, numseq, alnlen, gscores, 0.7, map, 'fasta', temp_res, "protein")
	Pal2Nal(temp_res, rawnuc, 'fasta', 'fasta', finalr_guidance70, 'fasta')
	shutil.copy(finalr_guidance70, '../'+alndir_nuc)
	shutil.copy(temp_res, '../'+alndir_aa+'/'+finalr_guidance70)
	
	#column masking
	mmod.maskColumns(refaln_file, numseq, alnlen, gscores, 0.7, map, 'fasta', temp_col, "protein")
	Pal2Nal(temp_col, rawnuc, 'fasta', 'fasta', finalc_guidance70, 'fasta')
	shutil.copy(finalc_guidance70, '../'+alndir_nuc)
	shutil.copy(temp_col, '../'+alndir_aa+'/'+finalc_guidance70)
	
		
		
	#######Weighted guidance
	
	##Mask at 90
	
	#residue masking
	mmod.maskResidues(refaln_file, numseq, alnlen, gweightscores, 0.9, map, 'fasta', temp_res, "protein")
	Pal2Nal(temp_res, rawnuc, 'fasta', 'fasta', finalr_gweights90, 'fasta')
	shutil.copy(finalr_gweights90, '../'+alndir_nuc)
	shutil.copy(temp_res, '../'+alndir_aa+'/'+finalr_gweights90)
	
	#column masking. 
	mmod.maskColumns(refaln_file, numseq, alnlen, gweightscores, 0.9, map, 'fasta', temp_col, "protein")
	Pal2Nal(temp_col, rawnuc, 'fasta', 'fasta', finalc_gweights90, 'fasta')
	shutil.copy(finalc_gweights90, '../'+alndir_nuc)
	shutil.copy(temp_col, '../'+alndir_aa+'/'+finalc_gweights90)
	
	##Mask at 70
	
	#residue masking
	mmod.maskResidues(refaln_file, numseq, alnlen, gweightscores, 0.7, map, 'fasta', temp_res, "protein")
	Pal2Nal(temp_res, rawnuc, 'fasta', 'fasta', finalr_gweights70, 'fasta')
	shutil.copy(finalr_gweights70, '../'+alndir_nuc)
	shutil.copy(temp_res, '../'+alndir_aa+'/'+finalr_gweights70)
	
	#column masking
	mmod.maskColumns(refaln_file, numseq, alnlen, gweightscores, 0.7, map, 'fasta', temp_col, "protein")
	Pal2Nal(temp_col, rawnuc, 'fasta', 'fasta', finalc_gweights70, 'fasta')
	shutil.copy(finalc_gweights70, '../'+alndir_nuc)
	shutil.copy(temp_col, '../'+alndir_aa+'/'+finalc_gweights70)


			
	#Also copy reference alignment to alndir_aa and alndir+nuc. NEED TO UNMAP IT FIRST!!!
	outref='refaln'+str(blah)+'.fasta'
	parsed = AlignIO.read(refaln_file, 'fasta')
	refMSA=MultipleSeqAlignment([])
	row=0
	for record in parsed:
		aln_record=SeqRecord(record.seq,id=str(map[row+1]), description='')
		refMSA.append(aln_record)
		row+=1
	outhandle=open(outref, 'w')
	outhandle.write(refMSA.format('fasta'))
	outhandle.close()
	shutil.copy(outref, '../'+alndir_aa)
	Pal2Nal(outref, rawnuc, 'fasta', 'fasta', '../'+alndir_nuc+'/'+outref, 'fasta')


	# Remove files in BootDir
	vomit=os.listdir('.')
	for file in vomit:
		os.remove(file)
	os.chdir('../')
	os.remove(raw)
	os.remove(refaln_file)
	os.remove(prealn_file)
	os.remove('prealn.dnd') ##this is clustal specific

##### NOW ALL THE ALIGNMENTS ARE COMPLETED. WE CAN BEGIN TO BUILD THE TREES!

prefix=['res70_guidance', 'res90_guidance', 'res70_gweights', 'res90_gweights', 'col70_guidance', 'col90_guidance', 'col70_gweights', 'col90_gweights', 'refaln']
for x in range(100):
	for case in prefix:
		alnfile=alndir_aa+'/'+case+str(x)+'.fasta'
		newfile=treedir+'/'+case+str(x)+'.tre'
		command2='/share/apps/fasttree-2.1.3/FastTreeMP -nosupport -gamma -mlacc 2 -slownni -spr 4 '+alnfile+' > '+newfile
		runit=subprocess.call(command2, shell=True)

		## Check that the tree was made, and then try to make it until it's made.
		size=os.path.getsize(newfile)
		while size==0:
			runit=subprocess.call(command2, shell=True)
			size=os.path.getsize(newfile)
		
