## Made on 6/13 to run fubar on reference alignments, true trees
## qsub file is trfu_rescol.py
## In the qsub file, do an array job 1-24:1 to cycle through all the basedir
## Also need to sed which aligner we're using, in the qsub file as well (since return directory also needs editing)

import re, os, sys, subprocess, fnmatch
import shutil

def GetFiles(ext, dirpath):
        files=[]
        dir=os.listdir(dirpath)
        for file in dir:
                if fnmatch.fnmatch(file, str(ext)+'*'):
                        files.append(file)
        return files
########

rundir='FUBARmaterials/'

## first entry of alndir empty since index 0
alndir=['', 'shortindel_64r', 'medindel_64r']
basedir=alndir[int(sys.argv[1])]

os.mkdir('seqs/')

#copy alignments
command1='cp /home/sjs3495/mafft_alntree/nucguided_'+basedir+'/refaln* seqs/'
call1=subprocess.call(command1, shell=True)


## copy truetree
command2='cp /home/sjs3495/current/truetrees/64r.tre '+rundir+'tree.tre'
call2=subprocess.call(command2, shell=True)
	
	
outdir='reffu_'+basedir
os.mkdir(outdir)
os.chdir(rundir)
for n in range(14,100):
	alnfile='../seqs/refaln'+str(n)+'.fasta'
	shutil.copy(alnfile, '.')
	shutil.move('refaln'+str(n)+'.fasta', 'temp.fasta')		
				
	cline='/home/sjs3495/bin/bin/HYPHYMP autoFUBAR.bf CPU=1'
	runit=subprocess.call(cline, shell=True)
	final_file='tree.tre.fubar.csv'
	shutil.move(final_file, '../'+outdir+'/refaln'+str(n)+'.fubar')
	
	#remove the plethora of fubar vomit
	vomit=os.listdir('.')
	for output in vomit:
		if fnmatch.fnmatch(output, "tree.tre.*"):
			os.remove(output)
	os.remove('temp.fasta')



















