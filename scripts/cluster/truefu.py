## 8/6/13


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

n=int(sys.argv[1])
if n==100:
	n=0
gene=sys.argv[2]
base=sys.argv[3]


command1='cp /home/sjs3495/rawdata_HA/'+gene+'/'+base+'/truealn_codon'+str(n)+'.fasta '+rundir
call1=subprocess.call(command1, shell=True)

## Copy over the true tree
find=re.search('seqs_(\w)', base)
if find:
	letter=find.group(1)
	command='cp /home/sjs3495/rawdata_HA/truetrees/'+letter+'_'+gene+'.tre '+rundir+'tree.tre'
	subprocess.call(command, shell=True)

	outdir='truefu_'+gene+'_'+base+'/'
	os.mkdir(outdir)
	os.chdir(rundir)
	shutil.move('truealn_codon'+str(n)+'.fasta', 'temp.fasta')		
					
	cline='/home/sjs3495/bin/bin/HYPHYMP autoFUBAR.bf CPU=1'
	runit=subprocess.call(cline, shell=True)
	final_file='tree.tre.fubar.csv'
	shutil.move(final_file, '../'+outdir+'/truealn'+str(n)+'.fubar')



















