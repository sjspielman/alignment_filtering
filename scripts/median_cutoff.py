import numpy as np
import re, csv

sweeps=csv.reader(open('fubar_sweeps.txt','r'))


count_old=0
temp=[]

cuts=[]

for row in sweeps:

	if row[3]=='or5':

		count=row[0] #count
		print count
		if int(count)==count_old: # if we're still in that count
			fp=float(row[2]) #fprate

			if (fp>0.0099 and fp<0.011):
				temp.append(float(row[1])) #cutoff
				print row[1]
		else:
			# get the last median and add to list
			cuts.append(np.median(temp))
			count_old=row[0]
				
			temp=[]
			fp=row[2] #fprate
			if (fp>0.0099 and fp<0.011):
				temp.append(row[1]) #cutoff
	else:
		break
		
print cuts
print np.median(cuts)
		

