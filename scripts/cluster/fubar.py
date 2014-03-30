## Rewritten 3/29/14

import os, sys, subprocess, fnmatch, shutil

rundir='FUBARmaterials/'

n=int(sys.argv[1]) - 1
gene=sys.argv[2]
alndir=sys.argv[3]+'/'
treedir=sys.argv[4]+'/'

outdir='fubar_'+gene+'/'
os.mkdir(outdir)


masks=['30_','50_', '70_', '90_']
algs=['Guidance_', 'BMweights_', 'PDweights_', 'GuidanceP_', 'BMweightsP_', 'PDweightsP_', 'refaln']


shutil.copy(treedir+"aatree"+str(n)+'.txt', rundir+'tree.tre')


os.chdir(rundir)
for alg in algs:
	for x in masks:
	
		if alg=='refaln':
			name='refaln'+str(n)+'.fasta'
		else:
			name=alg+x+str(n)+".fasta"	
		shutil.copy(alndir+name, ".")

		shutil.move(name, 'temp.fasta')

		cline='/home/sjs3495/bin/bin/HYPHYMP autoFUBAR.bf'
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




















