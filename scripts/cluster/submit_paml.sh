######## Submission script for paml jobs, via $QSUB only. Python script no longer needed.
## DO NOT RUN FLAT WITH PAML!!


## Starting conditions:
#GENE=or5
#MASK=50
#ALG=refaln

QSUB=/home/sjs3495/alignment_filtering/scripts/cluster/fubar.qsub

## Run reference alignment
qsub $QSUB

####### Run or5
sed -i 's/ALG=refaln/ALG=Guidance/g'
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
sed -i 's/ALG=PDweightsP/ALG=Guidance/g' $QSUB
sed -i 's/GENE=or5/GENE=rho/g' $QSUB

######## Run rho
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
sed -i 's/ALG=PDweightsP/ALG=Guidance/g' $QSUB
sed -i 's/GENE=rho/GENE=prk/g' $QSUB

######## Run prk
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


################### RUN GUIDANCE(P) FOR 30,70,90 MASKS ########################

## Reset algorithm to Guidance
sed -i 's/ALG=PDweightsP/ALG=Guidance/g' $QSUB

######### MASK 30
sed -i 's/MASK=50/MASK=30/g' $QSUB

## or5 
sed -i 's/GENE=prk/GENE=or5/g' $QSUB
qsub $QSUB
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=GuidanceP/ALG=Guidance/g' $QSUB

#rho
sed -i 's/GENE=or5/GENE=rho/g' $QSUB
qsub $QSUB
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=GuidanceP/ALG=Guidance/g' $QSUB

#prk
sed -i 's/GENE=rho/GENE=prk/g' $QSUB
qsub $QSUB
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=GuidanceP/ALG=Guidance/g' $QSUB


######### MASK 70
sed -i 's/MASK=30/MASK=70/g' $QSUB

## or5 
sed -i 's/GENE=prk/GENE=or5/g' $QSUB
qsub $QSUB
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=GuidanceP/ALG=Guidance/g' $QSUB

#rho
sed -i 's/GENE=or5/GENE=rho/g' $QSUB
qsub $QSUB
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=GuidanceP/ALG=Guidance/g' $QSUB

#prk
sed -i 's/GENE=rho/GENE=prk/g' $QSUB
qsub $QSUB
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=GuidanceP/ALG=Guidance/g' $QSUB


######### MASK 90
sed -i 's/MASK=70/MASK=90/g' $QSUB

## or5 
sed -i 's/GENE=prk/GENE=or5/g' $QSUB
qsub $QSUB
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=GuidanceP/ALG=Guidance/g' $QSUB

#rho
sed -i 's/GENE=or5/GENE=rho/g' $QSUB
qsub $QSUB
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=GuidanceP/ALG=Guidance/g' $QSUB

#prk
sed -i 's/GENE=rho/GENE=prk/g' $QSUB
qsub $QSUB
sed -i 's/ALG=Guidance/ALG=GuidanceP/g' $QSUB
qsub $QSUB
sed -i 's/ALG=GuidanceP/ALG=Guidance/g' $QSUB


########################### FULL RESET TO INITIAL CONDITIONS #############################

sed -i 's/MASK=90/MASK=50/g' $QSUB
sed -i 's/GENE=prk/GENE=or5/g' $QSUB
sed -i 's/ALG=Guidance/ALG=refaln/g' $QSUB


