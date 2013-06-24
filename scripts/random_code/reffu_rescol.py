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

aligner='clustal/'

rundir='FUBARmaterials/'

## first entry of alndir empty since index 0
alndir=['', 'shortindel_16p', 'medindel_16p', 'bigindel_16p', 'shortindel_16b', 'medindel_16b', 'bigindel_16b', 'shortindel_16r', 'medindel_16r', 'bigindel_16r', 'shortindel_16rc', 'medindel_16rc', 'bigindel_16rc', 'shortindel_64p', 'medindel_64p', 'bigindel_64p', 'shortindel_64b', 'medindel_64b', 'bigindel_64b', 'shortindel_64r', 'medindel_64r', 'bigindel_64r', 'shortindel_64rc', 'medindel_64rc', 'bigindel_64rc']
num=int(sys.argv[1])
basedir=alndir[num]
aligner=str(sys.argv[2])

os.mkdir('seqs/')
command1='cp /home/sjs3495/'+aligner+'_alntree/nucguided_'+basedir+'/refaln* seqs/'
call1=subprocess.call(command1, shell=True)


## Some obscene hardcoding for getting truetree
if num==1 or num==2 or num==3:
	command2='cp /home/sjs3495/current/truetrees/16p.tre '+rundir+'tree.tre'
elif num==4 or num==5 or num==6:
	command2='cp /home/sjs3495/current/truetrees/16b.tre '+rundir+'tree.tre'
elif num==7 or num==8 or num==9:
	command2='cp /home/sjs3495/current/truetrees/16r.tre '+rundir+'tree.tre'
elif num==10 or num==11 or num==12:
	command2='cp /home/sjs3495/current/truetrees/16rc.tre '+rundir+'tree.tre'
elif num==13 or num==14 or num==15:
	command2='cp /home/sjs3495/current/truetrees/64p.tre '+rundir+'tree.tre'
elif num==16 or num==17 or num==18:
	command2='cp /home/sjs3495/current/truetrees/64b.tre '+rundir+'tree.tre'
elif num==19 or num==20 or num==21:
	command2='cp /home/sjs3495/current/truetrees/64r.tre '+rundir+'tree.tre'
elif num==22 or num==23 or num==24:
	command2='cp /home/sjs3495/current/truetrees/64rc.tre '+rundir+'tree.tre'

call2=subprocess.call(command2, shell=True)
	
	
outdir='reffu_'+basedir
os.mkdir(outdir)
os.chdir(rundir)
for n in range(100):
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



















