## 8/6/13

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

n=int(sys.argv[1])
if n==100:
	n=0
gene=sys.argv[2]
base=sys.argv[3]

outdir='fubar_'+gene+'_'+base+'/'
os.mkdir(outdir)

alndir='cp /home/sjs3495/alntree/nucguided_linsi_'+gene+'_'+base+'/'
treedir='cp /home/sjs3495/alntree/aatrees_linsi_'+gene+'_'+base+'/aatree'

masks=['30_']#['50_', '70_', '90_']
algs=['guidance', 'BMweights', 'PDweights', 'guidance_p', 'BMweights_p', 'PDweights_p']#, 'refaln']


copy=treedir+str(n)+'.txt '+rundir+'tree.tre'
subprocess.call(copy, shell=True)

os.chdir(rundir)
for alg in algs:
	for x in masks:
	
		if alg=='refaln':
			name='refaln'+str(n)+'.fasta'
			copy=alndir+name+' .'
		else:
			name=alg+x+str(n)+".fasta"
			copy=alndir+name+' .'
		
		subprocess.call(copy, shell=True)

		shutil.move(name, 'temp.fasta')

		cline='/home/sjs3495/bin/bin/HYPHYMP autoFUBAR.bf CPU=1'
		runit=subprocess.call(cline, shell=True)
		final_file='tree.tre.fubar.csv'
		shutil.move(final_file, '../'+outdir+name+'.fubar')
	
		#remove the plethora of fubar vomit
		vomit=os.listdir('.')
		for output in vomit:
			if fnmatch.fnmatch(output, "tree.tre.*"):
				os.remove(output)
		os.remove('temp.fasta')
	
		if alg=='refaln': # just stop everything after refaln was processed
			assert 1==0




















