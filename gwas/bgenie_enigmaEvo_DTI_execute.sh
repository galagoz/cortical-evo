#----------------------------------------------------------------------
#!/bin/sh
#$ -N bgenie_enigmaEvo_DTI
#$ -q single15.q
#$ -S /bin/bash
#$ -e /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/log
#$ -o /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/log
#$ -m eas
#$ -M Barbara.Molz@mpi.nl
#----------------------------------------------------------------------

# this script runs BGENIE for different samples on global Freesurfer metrics, which are only controlled for the standard covariates
#written by B.Molz Nov 2021
chr=$1
sample=$2


##--------------------
date +"%D %R"

echo 'Running bgenie on chromosome' ${chr} 'for sample' ${sample} 

#------------------------------------------------------------------------------------------
# set up stuff e.g. output directory and directory where your bgen files are sitting
bgen_dir=/data/clusterfs/lag/users/barmol/enigma_evo/output/subset_dti
working_dir=/data/clusterfs/lag/users/barmol/enigma_evo/output/GWAS/dti
bgenie_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/DTI/GWAS

if [ ! -d ${working_dir} ] && echo "Directory DOES NOT exists - creating it ."; then
	mkdir -p ${working_dir}
fi
cd $working_dir

#------------------------------------------------------------------------------------------
# set up BGENIE input files
bgen_f=ukb43760_enigmaEvo_${sample}_chr${chr}.bgen
pheno_f=ukb43670_DTI_phenotype_bgenie.table
covar_f=ukb43670_DTI_confound_bgenie.table
snp_f=ukb43760_enigmaEvo_${sample}_snps2keep_chr${chr}.table
out_f=bgenie_ukb43760_${sample}_chr${chr}.out

#-----------------------------------------------------------------------------------------------------------------------------------------------
# run BGENIE on phenotype files without global FS correction on both global and regional values
if [ ! -f $working_dir/$out_f ]
then
	/data/workspaces/lag/shared_spaces/Resource_DB/bgenie/v1.2/bgenie --bgen $bgen_dir/$bgen_f \
			--thread 8 \
			--pvals \
			--include_rsids $bgenie_dir/$snp_f \
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