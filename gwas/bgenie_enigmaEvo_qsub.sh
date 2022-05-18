#!/bin/sh
#-----------------------------------------------------------------------------
# WRAPPER script: submit all BGENIE jobs to the cluster:
# Will run BGENIE on global FS metrics for different samples (Replication, European, different metrics and covariates)
# Written and adapted by B. Molz, Nov 2021
#------------------------------------------------------------------------------

for chrom in {1..21}; do
	for sample in replication european; do
	# lets run GLOBAL metrics for the respecitive sample (this has NO specific covariates)
		echo "RUN qsub bgenie_enigmaEvo_global_execute ${chrom} ukb43760_globalValues ${sample} "
		qsub -N bgenie_${sample}_global_${chrom} bgenie_enigmaEvo_global_execute.sh ${chrom} ${sample} 
	done
done

#run regional for replication 
for pheno in surface thickness; do
	for covar in noGlob withGlob; do 
		for chrom in {1..21}; do
			echo "RUN qsub bgenie_enigmaEvo_regional_replication_execute ${chrom} ukb43760_ regionalDK ${sample}"
			qsub -N bgenie_replication_${pheno}_${covar}_${chrom} bgenie_enigmaEvo_regional_replication_execute.sh ${chrom} $pheno $covar 
		done
	done
done


#run regional for european hemishpere 
for pheno in surface thickness; do
	for hemi in le re; do
		for covar in noGlob withGlob; do 
			for chrom in {1..21}; do
			# lets run regional metrics for the respecitive sample (this has two sets of covariates per sample, with global FS and without Global FS)
				echo "RUN qsub bgenie_enigmaEvo_regional_execute ${chrom} ukb43760_ regionalDK  $pheno $hemi $covar european"
				qsub -N bgenie_european_${pheno}_${hemi}_${covar}_${chrom} bgenie_enigmaEvo_regional_european_execute.sh ${chrom}  $pheno $hemi $covar
			done
		done
	done
done
