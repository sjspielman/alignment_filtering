## 6/24/13
#starts with shortindel_16p and clustal

# clustal
qsub paml.qsub
sed -i 's/shortindel_16p/medindel_16p/g' paml.qsub
qsub paml.qsub
sed -i 's/medindel_16p/bigindel_16p/g' paml.qsub
qsub paml.qsub
sed -i 's/bigindel_16p/shortindel_64p/g' paml.qsub
qsub paml.qsub
sed -i 's/shortindel_64p/medindel_64p/g' paml.qsub
qsub paml.qsub
sed -i 's/medindel_64p/bigindel_64p/g' paml.qsub
qsub paml.qsub
#reset basedir
sed -i 's/bigindel_64p/shortindel_16p/g' paml.qsub

# mafft
sed -i 's/clustal/mafft/g' paml.qsub
qsub paml.qsub
sed -i 's/shortindel_16p/medindel_16p/g' paml.qsub
qsub paml.qsub
sed -i 's/medindel_16p/bigindel_16p/g' paml.qsub
qsub paml.qsub
sed -i 's/bigindel_16p/shortindel_64p/g' paml.qsub
qsub paml.qsub
sed -i 's/shortindel_64p/medindel_64p/g' paml.qsub
qsub paml.qsub
sed -i 's/medindel_64p/bigindel_64p/g' paml.qsub
qsub paml.qsub
#reset basedir
sed -i 's/bigindel_64p/shortindel_16p/g' paml.qsub

# linsi
sed -i 's/mafft/linsi/g' paml.qsub
qsub paml.qsub
sed -i 's/shortindel_16p/medindel_16p/g' paml.qsub
qsub paml.qsub
sed -i 's/medindel_16p/bigindel_16p/g' paml.qsub
qsub paml.qsub
sed -i 's/bigindel_16p/shortindel_64p/g' paml.qsub
qsub paml.qsub
sed -i 's/shortindel_64p/medindel_64p/g' paml.qsub
qsub paml.qsub
sed -i 's/medindel_64p/bigindel_64p/g' paml.qsub
qsub paml.qsub

#reset basedir, aligner
sed -i 's/linsi/clustal/g' paml.qsub
sed -i 's/bigindel_64p/shortindel_16p/g' paml.qsub