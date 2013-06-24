## 6/13/13
## accompanies TAU_fubar.qsub and runTAU_fubar.py
## begins with aligner=clustal and basedir=shortindel_16p


qsub TAU_fubar.qsub #run p
sed -i 's/shortindel_16p/shortindel_16b/g' TAU_fubar.qsub
qsub TAU_fubar.qsub #run b
sed -i 's/shortindel_16b/shortindel_16r/g' TAU_fubar.qsub
qsub TAU_fubar.qsub #run r
sed -i 's/shortindel_16r/shortindel_16rc/g' TAU_fubar.qsub
qsub TAU_fubar.qsub #run rc

#switch aligner to mafft
sed -i 's/clustal/mafft/g' TAU_fubar.qsub
sed -i 's/shortindel_16rc/shortindel_16p/g' TAU_fubar.qsub
qsub TAU_fubar.qsub #run p
sed -i 's/shortindel_16p/shortindel_16b/g' TAU_fubar.qsub
qsub TAU_fubar.qsub #run b
sed -i 's/shortindel_16b/shortindel_16r/g' TAU_fubar.qsub
qsub TAU_fubar.qsub #run r
sed -i 's/shortindel_16r/shortindel_16rc/g' TAU_fubar.qsub
qsub TAU_fubar.qsub #run rc

# switch all back to original settings
sed -i 's/bigindel_64rc/shortindel_16p/g' TAU_fubar.qsub
sed -i 's/mafft/clustal/g' TAU_fubar.qsub
