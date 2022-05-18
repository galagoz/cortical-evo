#!/bin/sh
#----------------------------------------------------------------
# This scrip runs LDSC genetic correlations between derived summary statistics and UKBB derived summary statistics
#Info about LDSC: https://github.com/bulik/ldsc
#Info about summary stats from UKBB: https://open.win.ox.ac.uk/ukbiobank/big40/
# by Barbara Molz, Nov 2021
#----------------------------------------------------------------------------------------------------------
#set up the python environment needed on lux13 for LDSC
module purge
module load miniconda/3.2021.10 ldsc/v1.0.1
source activate ldsc

#set date 
DATE=`date +"%Y%m%d"`
#---------------------------------------------------------------------------------------------------
#set up our files: the general input path and the outputfolder
#and munged sumstats
#---------------------------------------------------------------------------------------------------
ukb_munged=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/DTI/LDSC/munged_UKBB
DTI_munged=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/DTI/LDSC/munged_enigmaEvo_DTI

#output gen cor
out_path_gencor=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/DTI/LDSC/gencor
# details for LDSC
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB
ldscores_dir=${resource_dir}/LDscores

#---------------------------------------------------------------------------------------------------
#set up stuff for gen cor
#----------------------------------------------------------------------------------------------------
mkdir -p ${out_path_gencor}
cd $ukb_munged
pwd
for m_sumstats in ${DTI_munged}/*gz; do
	filename=${m_sumstats##*/}
	pheno=${filename##*europeanDTI_}
	roi=${pheno%_allChr*gz}
	echo 'curren roi' $roi
	e3_file=$(find -maxdepth 1 -type f -iname "*${roi}*.gz" -printf "%f\n")
	echo 'curren ukbb:' $e3_file
	if [ ! -f ${out_path_gencor}/gencor_ukb43760_DTI_${roi}.log ]; then
		echo 'run LDSC to calculate genetic correlation between E3 and enigma Evo' ${e3_file} ${roi}
		ldsc.py \
		--rg ${e3_file},${m_sumstats} \
		--ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
		--w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
		--out ${out_path_gencor}/gencor_ukb43760_DTI_${roi}
	fi
done


# ---------------------------------------------------------------------------------------------------
#extract gen cor resuls 
#----------------------------------------------------------------------------------------------------
#little wrapper to extract the summary values from the LDSC log after running gen cor
echo "trait rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se"> $out_path_gencor/gencor_ukbBig40_DTI.txt
for i in $out_path_gencor/*log; do
  res=$(cat $i | grep '^m.*sumstat' | cut -d" " -f4-)
  pheno=$(basename "$i" | cut -d. -f1)
  echo $pheno $res >> $out_path_gencor/gencor_ukbBig40_DTI.txt
done 