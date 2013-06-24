### linsi runs
qsub custom_fubar_rescol.qsub
sed -i "s/col90_guidance/res70_gweights/g" custom_fubar_rescol.py
sed -i "s/x=19/x=96/g" custom_fubar_rescol.py
qsub custom_fubar_rescol.qsub
sed -i "s/res70_gweights/refaln/g" custom_fubar_rescol.py
sed -i "s/x=96/x=81/g" custom_fubar_rescol.py
qsub custom_fubar_rescol.qsub
sed -i "s/bigindel_64p/shortindel_64rc/g" custom_fubar_rescol.py
sed -i "s/x=81/x=44/g" custom_fubar_rescol.py
qsub custom_fubar_rescol.qsub
sed -i "s/shortindel_64rc/medindel_64rc/g" custom_fubar_rescol.py
sed -i "s/x=44/x=76/g" custom_fubar_rescol.py
qsub custom_fubar_rescol.qsub
sed -i "s/medindel_64rc/bigindel_64rc/g" custom_fubar_rescol.py
sed -i "s/refaln/col70_gweights/g" custom_fubar_rescol.py
sed -i "s/x=76/x=75/g" custom_fubar_rescol.py
qsub custom_fubar_rescol.qsub

## mafft run
sed -i "s/linsi/mafft/g" custom_fubar_rescol.py
sed -i "s/linsi/mafft/g" custom_fubar_rescol.qsub
sed -i "s/bigindel_64rc/medindel_64rc/g" custom_fubar_rescol.py
sed -i "s/col70_gweights/col70_guidance/g" custom_fubar_rescol.py
sed -i "s/x=75/x=31/g" custom_fubar_rescol.py
qsub custom_fubar_rescol.qsub

##clustal runs
sed -i "s/mafft/clustal/g" custom_fubar_rescol.py
sed -i "s/mafft/clustal/g" custom_fubar_rescol.qsub
sed -i "s/medindel_64rc/medindel_64b/g" custom_fubar_rescol.py
sed -i "s/col70_guidance/col90_gweights/g" custom_fubar_rescol.py
sed -i "s/x=31/x=65/g" custom_fubar_rescol.py
qsub custom_fubar_rescol.qsub


sed -i "s/medindel_64b/medindel_64p/g" custom_fubar_rescol.py
sed -i "s/x=65/x=82/g" custom_fubar_rescol.py
qsub custom_fubar_rescol.qsub

sed -i "s/medindel_64p/bigindel_64b/g" custom_fubar_rescol.py
sed -i "s/col90_gweights/res90_gweights/g" custom_fubar_rescol.py
sed -i "s/x=82/x=88/g" custom_fubar_rescol.py
qsub custom_fubar_rescol.qsub

sed -i "s/bigindel_64b/shortindel_64rc/g" custom_fubar_rescol.py
sed -i "s/x=88/x=95/g" custom_fubar_rescol.py
qsub custom_fubar_rescol.qsub

sed -i "s/shortindel_64rc/bigindel_64rc/g" custom_fubar_rescol.py
sed -i "s/res90_gweights/res70_gweights/g" custom_fubar_rescol.py
sed -i "s/x=95/x=22/g" custom_fubar_rescol.py
qsub custom_fubar_rescol.qsub

## RESET!
sed -i "s/clustal/linsi/g" custom_fubar_rescol.py
sed -i "s/clustal/linsi/g" custom_fubar_rescol.qsub
sed -i "s/bigindel_64rc/bigindel_64p/g" custom_fubar_rescol.py
sed -i "s/res70_gweights/col90_guidance/g" custom_fubar_rescol.py
sed -i "s/x=22/x=19/g" custom_fubar_rescol.py
