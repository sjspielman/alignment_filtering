######## Submission script for alntree jobs, via $QSUB only. Simply needs to cycle across genes #########

QSUB=/home/sjs3495/alignment_filtering/scripts/cluster/alntree.qsub
# starts with GENE=or5

## or5
qsub $QSUB

## rho
sed -i 's/GENE=or5/GENE=rho/g' $QSUB
qsub $QSUB

## prk
sed -i 's/GENE=rho/GENE=prk/g' $QSUB
qsub $QSUB

## flat
sed -i 's/GENE=prk/GENE=flat/g' $QSUB
qsub $QSUB

## RESET TO or5
sed -i 's/GENE=flat/GENE=or5/g' $QSUB
