### 6/14/13
### Accompanies apply_optimal.py and optimal.qsub
### Starts w/ clustal

qsub optimal.qsub
sed -i "s/clustal/mafft/g" optimal.qsub
qsub optimal.qsub
sed -i "s/mafft/linsi/g" optimal.qsub
qsub optimal.qsub
sed -i "s/linsi/clustal/g" optimal.qsub
