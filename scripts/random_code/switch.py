

bases=['shortindel_16p', 'medindel_16p', 'bigindel_16p', 'shortindel_16b', 'medindel_16b', 'bigindel_16b', 'shortindel_16r', 'medindel_16r', 'bigindel_16r', 'shortindel_16rc', 'medindel_16rc', 'bigindel_16rc', 'shortindel_64p', 'medindel_64p', 'bigindel_64p', 'shortindel_64b', 'medindel_64b', 'bigindel_64b', 'shortindel_64r', 'medindel_64r', 'bigindel_64r', 'shortindel_64rc', 'medindel_64rc', 'bigindel_64rc']

for n in range(1,len(bases)):
	command="sed -i 's/"+bases[n-1]+"/"+bases[n]+"/g' truefubar.qsub"
	print command
	print "qsub truefubar.qsub"