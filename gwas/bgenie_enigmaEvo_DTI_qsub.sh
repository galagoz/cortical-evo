#!/bin/sh

#WRAPPER script: submit script for DTI sample to cluster
# Will run BGENIE on global FS metrics for different samples (here European all / LERE/ replication)

for chrom in {1..20}; do
	echo "RUN bgenie_enigmaEvo_global_execute.sh ${chrom} ukb43760_enigmaEvo DTI"
	qsub  bgenie_enigmaEvo_DTI_execute.sh ${chrom} europeanDTI
done
