## 6/24/13. 

import re, os, sys, subprocess, fnmatch
import shutil
from Bio import AlignIO

def GetFiles(ext, dirpath):
        files=[]
        dir=os.listdir(dirpath)
        for file in dir:
                if fnmatch.fnmatch(file, str(ext)+'.+'):
                        files.append(file)
        return files
########


## Deal with the base in the qsub file.

filedict=['', 'res70_guidance', 'res90_guidance', 'res70_gweights', 'res90_gweights', 'refaln']

prefix=filedict[int(sys.argv[1])]
basedir=str(sys.argv[2]) ## eg, shortindel_16p, medindel_16p, bigindel_16p, and 64 of those.
aligner=str(sys.argv[3]) ## eg, clustal

alndir='alns/'
treedir='trees/'
os.mkdir(alndir)
os.mkdir(treedir)


# Copy files over
for n in range(50):
	file1='/home/sjs3495/current/all_results/'+aligner+'_alntree/nucguided_'+basedir+'/'+prefix+str(n)+'.fasta'
	command1='cp '+file1+' '+alndir
	print command1
	run1=subprocess.call(command1, shell=True)
	file2='/home/sjs3495/current/all_results/'+aligner+'_alntree/aatrees_'+basedir+'/'+prefix+str(n)+'.tre'
	command2='cp '+file2+' '+treedir
	print command2
	run2=subprocess.call(command2, shell=True)	
	
outdir='paml_'+basedir
os.mkdir(outdir)

for n in range(50):
	name=prefix+str(n)
	
	# Convert alignment to phylip
	alnfile=alndir+name+'.fasta'
	print alnfile
	out_alnfile=name+'.phy'
	parsed=(AlignIO.read(alnfile, 'fasta'))
	outhandle=open(out_alnfile, 'w')
	AlignIO.write(parsed, outhandle, 'phylip')
	outhandle.close()
	format = "sed -i '1s/$/ I/' "+out_alnfile
	runFormat=subprocess.call(format, shell=True)

	# Rename files to match codeml.ctl
	shutil.move(out_alnfile, 'temp.phy')
	treefile=treedir+name+'.tre'
	shutil.move(treefile, 'tree.tre')				

			
	cline='/home/sjs3495/bin/codeml codeml.ctl'
	runit=subprocess.call(cline, shell=True)
	shutil.move('rst', outdir+'/'+name+'.paml')

	os.remove('tree.tre')
	os.remove('temp.phy')

	
	
	
	
	
	
	
	
	
	























