#!/bin/sh

#WRAPPER script to submit subsetting jobs to the cluster

for sample in replication european; do 
	for chrom in {1..20} ; do
		echo "RUN qsub subsetting_imaging40k_subset_and_snpstats.sh -s subsetting_imaging40k_subset_and_snpstats_config_${sample}_ukb43760.txt -c ${chrom}"	
		qsub  subsetting_imaging40k_subset_and_snpstats.sh -s subsetting_imaging40k_subset_and_snpstats_config_${sample}_ukb43760.txt -c ${chrom}
	done
done
