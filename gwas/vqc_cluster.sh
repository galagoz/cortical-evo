#!/bin/sh

#WRAPPEr script to submit variant qc to the grid for different samples 

for sample in european replication; do 
	for chrom in {1..20} X; do
		echo "Run variant QC for chrom ${chrom} on sample ${sample}"
		qsub -N variantQC_${sample}_${chrom} vqc_execute.sh $chrom $sample
	done
done
