## 6/13/13. Goes with qsub file fubar_rescol.qsub. MAKE SURE TO sed ALIGNER AND BASEDIR FOR RUNS!!

import re, os, sys, subprocess, fnmatch
import shutil

def GetFiles(ext, dirpath):
        files=[]
        dir=os.listdir(dirpath)
        for file in dir:
                if fnmatch.fnmatch(file, str(ext)+'.+'):
                        files.append(file)
        return files
########

rundir='FUBARmaterials/'

basedir='bigindel_64p'
aligner='linsi'
prefix='col90_guidance'
x=19     ## range needed to run.

alndir='alns/'
treedir='trees/'
os.mkdir(alndir)
os.mkdir(treedir)

for n in range(x,100):
	file1='/home/sjs3495/'+aligner+'_alntree/nucguided_'+basedir+'/'+prefix+str(n)+'.fasta'
	command1='cp '+file1+' '+alndir
	print command1
	run1=subprocess.call(command1, shell=True)
	file2='/home/sjs3495/'+aligner+'_alntree/aatrees_'+basedir+'/'+prefix+str(n)+'.tre'
	command2='cp '+file2+' '+treedir
	print command2
	run2=subprocess.call(command2, shell=True)	

	
outdir='fubar_'+basedir
os.mkdir(outdir)

os.chdir(rundir)
for n in range(x,100):
	name=prefix+str(n)
	treefile='../'+treedir+name+'.tre'
	alnfile='../'+alndir+name+'.fasta'
	
	shutil.copy(alnfile, '.')
	shutil.move(name+'.fasta', 'temp.fasta')
	shutil.copy(treefile, '.')
	shutil.move(name+'.tre', 'tree.tre')				
			
	cline='/home/sjs3495/bin/bin/HYPHYMP autoFUBAR.bf CPU=1'
	runit=subprocess.call(cline, shell=True)
	final_file='tree.tre.fubar.csv'
	shutil.move(final_file, '../'+outdir+'/'+name+'.fubar')
	
	#remove the plethora of fubar vomit
	vomit=os.listdir('.')
	for output in vomit:
		if fnmatch.fnmatch(output, "tree.tre.*"):
			os.remove(output)
	os.remove('tree.tre')
	os.remove('temp.fasta')























