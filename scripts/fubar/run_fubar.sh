## 6/13/13
## accompanies fubar_rescol.qsub and runfubar_rescol.py
## begins with aligner=clustal and basedir=shortindel_16p


qsub fubar_rescol.qsub
sed -i 's/shortindel_16p/medindel_16p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_16p/bigindel_16p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_16p/shortindel_16b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_16b/medindel_16b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_16b/bigindel_16b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_16b/shortindel_16r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_16r/medindel_16r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_16r/bigindel_16r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_16r/shortindel_16rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_16rc/medindel_16rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_16rc/bigindel_16rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_16rc/shortindel_64p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_64p/medindel_64p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_64p/bigindel_64p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_64p/shortindel_64b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_64b/medindel_64b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_64b/bigindel_64b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_64b/shortindel_64r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_64r/medindel_64r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_64r/bigindel_64r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_64r/shortindel_64rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_64rc/medindel_64rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_64rc/bigindel_64rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_64rc/shortindel_16p/g' fubar_rescol.qsub

#switch aligner to mafft
sed -i 's/clustal/mafft/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_16p/medindel_16p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_16p/bigindel_16p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_16p/shortindel_16b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_16b/medindel_16b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_16b/bigindel_16b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_16b/shortindel_16r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_16r/medindel_16r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_16r/bigindel_16r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_16r/shortindel_16rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_16rc/medindel_16rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_16rc/bigindel_16rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_16rc/shortindel_64p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_64p/medindel_64p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_64p/bigindel_64p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_64p/shortindel_64b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_64b/medindel_64b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_64b/bigindel_64b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_64b/shortindel_64r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_64r/medindel_64r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_64r/bigindel_64r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_64r/shortindel_64rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_64rc/medindel_64rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_64rc/bigindel_64rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_64rc/shortindel_16p/g' fubar_rescol.qsub

#switch aligner to linsi
sed -i 's/mafft/linsi/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_16p/medindel_16p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_16p/bigindel_16p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_16p/shortindel_16b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_16b/medindel_16b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_16b/bigindel_16b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_16b/shortindel_16r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_16r/medindel_16r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_16r/bigindel_16r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_16r/shortindel_16rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_16rc/medindel_16rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_16rc/bigindel_16rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_16rc/shortindel_64p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_64p/medindel_64p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_64p/bigindel_64p/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_64p/shortindel_64b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_64b/medindel_64b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_64b/bigindel_64b/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_64b/shortindel_64r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_64r/medindel_64r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_64r/bigindel_64r/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/bigindel_64r/shortindel_64rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/shortindel_64rc/medindel_64rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub
sed -i 's/medindel_64rc/bigindel_64rc/g' fubar_rescol.qsub
qsub fubar_rescol.qsub

# switch all back to original settings
sed -i 's/bigindel_64rc/shortindel_16p/g' fubar_rescol.qsub
sed -i 's/linsi/clustal/g' fubar_rescol.qsub
