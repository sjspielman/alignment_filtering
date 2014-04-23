######## Submission script for fubar jobs, via $QSUB only. Python script no longer needed.

## Starting conditions:
#GENE=or5
#MASK=50
#ALG=refaln

QSUB=/home/sjs3495/alignment_filtering/scripts/cluster/fubar.qsub

## Run reference alignment
qsub $QSUB

####### Run or5
sed -i 's/ALG=refaln/ALG=Guidance/g' $QSUB
qsub $QSUB
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=GuidanceP/ALG=BMweights/g' $QSUB
qsub $QSUB
sed -i 's/ALG=BMweights/ALG=BMweightsP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=BMweightsP/ALG=PDweights/g' $QSUB
qsub $QSUB
sed -i 's/ALG=PDweights/ALG=PDweightsP/g' $QSUB
qsub $QSUB

######## Reset algorithm and change gene
sed -i 's/ALG=PDweightsP/ALG=refaln/g' $QSUB
sed -i 's/GENE=or5/GENE=rho/g' $QSUB

######## Run rho
qsub $QSUB
sed -i 's/ALG=refaln/ALG=Guidance/g' $QSUB
qsub $QSUB
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=GuidanceP/ALG=BMweights/g' $QSUB
qsub $QSUB
sed -i 's/ALG=BMweights/ALG=BMweightsP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=BMweightsP/ALG=PDweights/g' $QSUB
qsub $QSUB
sed -i 's/ALG=PDweights/ALG=PDweightsP/g' $QSUB
qsub $QSUB

######## Reset algorithm and change gene
sed -i 's/ALG=PDweightsP/ALG=refaln/g' $QSUB
sed -i 's/GENE=rho/GENE=prk/g' $QSUB

######## Run prk
qsub $QSUB
sed -i 's/ALG=refaln/ALG=Guidance/g' $QSUB
qsub $QSUB
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=GuidanceP/ALG=BMweights/g' $QSUB
qsub $QSUB
sed -i 's/ALG=BMweights/ALG=BMweightsP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=BMweightsP/ALG=PDweights/g' $QSUB
qsub $QSUB
sed -i 's/ALG=PDweights/ALG=PDweightsP/g' $QSUB
qsub $QSUB

######## Reset algorithm and change gene
sed -i 's/ALG=PDweightsP/ALG=refaln/g' $QSUB
sed -i 's/GENE=prk/GENE=flat/g' $QSUB

######## Run flat
qsub $QSUB
sed -i 's/ALG=refaln/ALG=Guidance/g' $QSUB
qsub $QSUB
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=GuidanceP/ALG=BMweights/g' $QSUB
qsub $QSUB
sed -i 's/ALG=BMweights/ALG=BMweightsP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=BMweightsP/ALG=PDweights/g' $QSUB
qsub $QSUB
sed -i 's/ALG=PDweights/ALG=PDweightsP/g' $QSUB
qsub $QSUB



########################### FULL RESET TO INITIAL CONDITIONS #############################

sed -i 's/GENE=flat/GENE=or5/g' $QSUB
sed -i 's/ALG=PDweightsP/ALG=refaln/g' $QSUB


