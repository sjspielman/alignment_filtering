######## Submission script for all truealn jobs, both fubar and paml. ##########

## Starting conditions:
#DATASET=HA
#GENE=or5

QSUB_FU=/home/sjs3495/alignment_filtering/scripts/cluster/truefu.qsub
QSUB_PAML=/home/sjs3495/alignment_filtering/scripts/cluster/truepaml.qsub

######## Run all genes for HA, then all genes for GP41 ########

###################################### HA ################################################
# HA, or5
qsub $QSUB_FU
qsub $QSUB_PAML

# HA, rho
sed -i 's/GENE=or5/GENE=rho/g' $QSUB_FU
sed -i 's/GENE=or5/GENE=rho/g' $QSUB_PAML
qsub $QSUB_FU
qsub $QSUB_PAML

# HA, prk
sed -i 's/GENE=rho/GENE=prk/g' $QSUB_FU
sed -i 's/GENE=rho/GENE=prk/g' $QSUB_PAML
qsub $QSUB_FU
qsub $QSUB_PAML

# HA, flat. FUBAR ONLY.
sed -i 's/GENE=prk/GENE=flat/g' $QSUB_FU
qsub $QSUB_FU
##########################################################################################


## RESET TO OR5 AND SWAP DATASET TO GP41
sed -i 's/GENE=flat/GENE=or5/g' $QSUB_FU
sed -i 's/GENE=prk/GENE=or5/g' $QSUB_PAML
sed -i 's/DATASET=HA/DATASET=GP41/g' $QSUB_FU
sed -i 's/DATASET=HA/DATASET=GP41/g' $QSUB_PAML

##################################### GP41 ###############################################

# GP41, or5
qsub $QSUB_FU
qsub $QSUB_PAML

# GP41, rho
sed -i 's/GENE=or5/GENE=rho/g' $QSUB_FU
sed -i 's/GENE=or5/GENE=rho/g' $QSUB_PAML
qsub $QSUB_FU
qsub $QSUB_PAML

# GP41, prk
sed -i 's/GENE=rho/GENE=prk/g' $QSUB_FU
sed -i 's/GENE=rho/GENE=prk/g' $QSUB_PAML
qsub $QSUB_FU
qsub $QSUB_PAML

# GP41, flat. FUBAR ONLY.
sed -i 's/GENE=prk/GENE=flat/g' $QSUB_FU
qsub $QSUB_FU
##########################################################################################


#### RESET TO STARTING CONDITIONS
sed -i 's/GENE=flat/GENE=or5/g' $QSUB_FU
sed -i 's/GENE=prk/GENE=or5/g' $QSUB_PAML
sed -i 's/DATASET=GP41/DATASET=HA/g' $QSUB_FU
sed -i 's/DATASET=GP41/DATASET=HA/g' $QSUB_PAML

