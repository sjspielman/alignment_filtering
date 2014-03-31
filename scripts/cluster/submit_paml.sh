######## Submission script for paml jobs, via paml.qsub only. Python script no longer needed.

## For now, we are just starting out with rho neutral. Just run all algs with 50 mask and rho fixed.

## Starting conditions:
#GENE=rho
#MASK=50
#ALG=Guidance

qsub paml.qsub
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' paml.qsub
qsub paml.qsub
sed -i 's/ALG=GuidanceP/ALG=BMweights/g' paml.qsub
qsub paml.qsub
sed -i 's/ALG=BMweights/ALG=BMweightsP/g' paml.qsub
qsub paml.qsub
sed -i 's/ALG=BMweightsP/ALG=PDweights/g' paml.qsub
qsub paml.qsub
sed -i 's/ALG=PDweights/ALG=PDweightsP/g' paml.qsub
qsub paml.qsub
sed -i 's/ALG=PDweightsP/ALG=refaln/g' paml.qsub
qsub paml.qsub

sed -i 's/ALG=refaln/ALG=Guidance/g' paml.qsub