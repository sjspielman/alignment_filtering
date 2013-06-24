## 6/13/13.
## accompanies trfu_rescol.qsub, truefu_rescol.py, reffu_rescol.py
## starts with aligner=clustal and script=truefu_rescol.py

## Run all the true, true
qsub trfu_rescol.qsub
sed -i 's/truefu/reffu/g' trfu_rescol.qsub
qsub trfu_rescol.qsub
sed -i 's/reffu/truefu/g' trfu_rescol.qsub



sed -i 's/clustal/mafft/g' trfu_rescol.qsub
qsub trfu_rescol.qsub
sed -i 's/mafft/linsi/g' trfu_rescol.qsub
qsub trfu_rescol.qsub
sed -i 's/linsi/clustal/g' trfu_rescol.qsub

## Run all the ref, true
sed -i 's/truefu/reffu/g' trfu_rescol.qsub
qsub trfu_rescol.qsub
sed -i 's/clustal/mafft/g' trfu_rescol.qsub
qsub trfu_rescol.qsub
sed -i 's/mafft/linsi/g' trfu_rescol.qsub
qsub trfu_rescol.qsub

## Back to original settings
sed -i 's/linsi/clustal/g' trfu_rescol.qsub
sed -i 's/linsi/clustal/g' trfu_rescol.qsub
