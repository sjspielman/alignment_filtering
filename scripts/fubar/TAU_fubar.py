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

############# qsub file settings.
aaprefixes=['', 'res_aa70', 'res_aa90', 'col_aa70', 'col_aa90']
nucprefixes=['', 'res_nuc70', 'res_nuc90', 'col_nuc70', 'col_nuc90']
plainprefixes=['', 'res70', 'res90', 'col70', 'col90']
aaprefix=aaprefixes[int(sys.argv[1])]
nucprefix=nucprefixes[int(sys.argv[1])]
plainprefix=plainprefixes[int(sys.argv[1])]

basedir=str(sys.argv[2]) ## eg, shortindel_16b
aligner=str(sys.argv[3]) ## eg, clustal
#############

alndir='alns/'
treedir='trees/'
os.mkdir(alndir)
os.mkdir(treedir)

for n in range(100):
	file1='/home/sjs3495/TAU_alntree/'+aligner+'/TAU_'+basedir+'/'+nucprefix+str(n)+'.fasta'
	command1='cp '+file1+' '+alndir
	print command1
	run1=subprocess.call(command1, shell=True)
	file1='/home/sjs3495/TAU_alntree/'+aligner+'/aatrees_TAU_'+basedir+'/'+aaprefix+str(n)+'.fasta'
	command2='cp '+file2+' '+treedir
	print command2
	run2=subprocess.call(command2, shell=True)	

	
outdir='fubar_'+basedir
os.mkdir(outdir)

os.chdir(rundir)
for n in range(100):
	name=plainprefix+'_'+str(n)
	treefile='../'+treedir+'/'+name+'.tre'
	alnfile='../'+alndir+'/'+name+'.fasta'
	
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























