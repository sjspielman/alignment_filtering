### Script to convert all the truealn's to be correct, godfuckingdamnit
import re, os, sys, subprocess, fnmatch
import shutil
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

def GetFiles(ext, dirpath):
	files=[]
	dir=os.listdir(dirpath)
	for file in dir:
		if fnmatch.fnmatch(file, '*.'+str(ext)):
			files.append(file)
	return files
########

##################
def parseTrue(line):
	find=re.search('^(t\d+)\s+([\w-]+)', line)
	if find:
		id=find.group(1)
		seq=find.group(2)
	return (id,seq)
##################

rundict={1:'med128b', 2:'med16p', 3:'med64r', 4:'short16b', 5:'short64p', 6:'tiny16b', 7:'tiny64p', 8:'med128p', 9:'med16r', 10:'short128b', 11:'short16p', 12:'short64r', 13:'tiny16p', 14:'tiny64r', 15:'med128r', 16:'med64b', 17:'short128p', 18:'short16r', 19:'tiny128p', 20:'tiny16r', 21:'med16b', 22:'med64p', 23:'short128r', 24:'short64b', 25:'tiny128r', 26:'tiny64b', 27:'tiny128b'}
print "hello"
for n in range(1,28):
	basedir='guided_'+rundict[n]
	print basedir
	findnum=re.search('[medshortiny]+(\d+)[prb]', basedir)
	if findnum:
		lenaln=int(findnum.group(1))

	
	#next two lines in case of cluster.
	#command1='cp -r /home/sjs3495/current/guided_4.7/'+alndir+' .'
	#run1=subprocess.call(command1, shell=True)


	files=GetFiles('phy',basedir)
	for file in files:
		find=re.search('(truealn\d+)\.phy', file)
		if find:
			print file
			outfile=basedir+'/'+find.group(1)+'.phy'
			inhandle=open(basedir+'/'+file, 'r')
			lines=inhandle.readlines()
			inhandle.close()
			
			##skip the first line since it's all phylip-y
			newMSA=MultipleSeqAlignment([])
			for i in range(1, lenaln+1):
				(id,seq)=parseTrue(lines[i])
				record=SeqRecord(Seq(str(seq), generic_dna), id=id, description='')
				newMSA.append(record)
			outhandle=open(outfile, 'w')
			umm=AlignIO.write(newMSA, outhandle, 'phylip')
			outhandle.close()