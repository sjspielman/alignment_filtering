from Bio import SeqIO

class Map:
	def __init__(self):
		'''initialization function'''
		return
		
	def ids2int(self, unaligned_aa, format, prealn_file):	
		print "Building map"		
		idmap={}		
		n=0
		
		infile=open(unaligned_aa, 'r')
		parsed = list(SeqIO.parse(infile, str(format)))
		infile.close()		
		
		out=open(prealn_file, 'w')
		
		for record in parsed:	
			n+=1
			idmap[int(n)] = str(record.id)
			seq=str(record.seq)
			out.write('>'+str(n)+'\n'+seq+'\n')
		out.close()
		return (idmap)
	


	