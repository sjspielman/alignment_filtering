## 8/6/13

import re, os, sys, subprocess, fnmatch
import shutil

n=int(sys.argv[1])
if n==100:
	n=0
gene=sys.argv[2]
base=sys.argv[3]


outdir='truepaml_'+gene+'_'+base+'/'
os.mkdir(outdir)

command1='cp /home/sjs3495/rawdata_HA/'+gene+'/'+base+'/truealn_codon'+str(n)+'.fasta temp.fasta'
call1=subprocess.call(command1, shell=True)

## Copy over the true tree
find=re.search('seqs_(\w)', base)
if find:
	letter=find.group(1)
	command='cp /home/sjs3495/rawdata_HA/truetrees/'+letter+'_'+gene+'.tre tree.tre'
	subprocess.call(command, shell=True)
					
	cline='/home/sjs3495/bin/codeml codeml.ctl'
	subprocess.call(cline, shell=True)
	final_file='rst'
	shutil.move(final_file, outdir+'truealn'+str(n)+'.rst')
















