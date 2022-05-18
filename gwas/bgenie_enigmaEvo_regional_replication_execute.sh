#----------------------------------------------------------------------
#!/bin/sh
#$ -q single15.q
#$ -S /bin/bash
#$ -e /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/log
#$ -o /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/log
#$ -m eas
#$ -M Barbara.Molz@mpi.nl
#----------------------------------------------------------------------

# this script runs BGENIE on the replication sample for regional Freesurfer metrics, which are controlled for two sets of covariates
# written and adapted by B. Molz, 2021
chr=$1
pheno=$2
covar=$3

##--------------------
date +"%D %R"

echo 'Running bgenie for replication on' ${pheno}' on chromosome' ${chr} 'with confounds' ${covar}

#------------------------------------------------------------------------------------------
# set up stuff e.g. output directory and directory where your bgen files are sitting
bgen_dir=/data/clusterfs/lag/users/barmol/enigma_evo/output/subset_replication
working_dir=/data/clusterfs/lag/users/barmol/enigma_evo/output/GWAS/replication
bgenie_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/replication/bgenieInput
cd $working_dir
 
#------------------------------------------------------------------------------------------
# set up BGENIE input files
bgen_f=ukb43760_enigmaEvo_replication_chr${chr}.bgen
pheno_f=ukb43760_regional_${pheno}_replication.table
snp_f=ukb43760_enigmaEvo_replication_snps2keep_chr${chr}.table
out_f=bgenie_ukb43760_regional_${pheno}_${covar}_chr${chr}_replication.out

if [[ ${covar} == "withGlob" ]]; then
	covar_f=ukb43760_covariates_${pheno}_replication.table
else
	covar_f=ukb43760_covariates_noGlob_replication.table
fi


#-----------------------------------------------------------------------------------------------------------------------------------------------
# run BGENIE on phenotype files without global FS correction on both global and regional values
if [ ! -f $working_dir/$out_f ];then
	/data/workspaces/lag/shared_spaces/Resource_DB/bgenie/v1.2/bgenie --bgen $bgen_dir/$bgen_f \
	--thread 8 \
	--pvals \
	--include_rsids $bgenie_dir/snps2keep/$snp_f \
	--pheno $bgenie_dir/$pheno_f \
	--covar $bgenie_dir/$covar_f \
	--out $working_dir/$out_f
else
	echo 'OUTPUT FILE ALREADY EXISTST, SKIP'
fi

date +"%D %R" 
#####
##
#