#----------------------------------------------------------------------
#!/bin/sh
#$ -q single15.q
#$ -S /bin/bash
#$ -e /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/log
#$ -o /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/log
#$ -m eas
#$ -M Barbara.Molz@mpi.nl
#----------------------------------------------------------------------

# this script runs BGENIE for different samples on global Freesurfer metrics, which are only controlled for the standard covariates
chr=$1
sample=$2

##--------------------
date +"%D %R"

echo 'Running bgenie on global Freesurfer metrics on chromosome' ${chr} 'for sample' ${sample} 

#------------------------------------------------------------------------------------------
# set up stuff e.g. output directory and directory where your bgen files are sitting
bgen_dir=/data/clusterfs/lag/users/barmol/enigma_evo/output/subset_${sample}
working_dir=/data/clusterfs/lag/users/barmol/enigma_evo/output/GWAS/${sample}
bgenie_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/$sample/bgenieInput
cd $working_dir

#------------------------------------------------------------------------------------------
# set up BGENIE input files
bgen_f=ukb43760_enigmaEvo_${sample}_chr${chr}.bgen
pheno_f=ukb43760_global_${sample}.table
covar_f=ukb43760_covariates_noGlob_${sample}.table
snp_f=ukb43760_enigmaEvo_${sample}_snps2keep_chr${chr}.table
out_f=bgenie_ukb43760_global_chr${chr}_${sample}.out

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