#!/bin/sh
#----------------------------------------------------------------# This scrip runs LDSC genetic correlations between derived DIT summary statistics to create correlation matrix for PhenoSPD
#Info about LDSC: https://github.com/bulik/ldsc
# by Barbara Molz, Nov 2021
#----------------------------------------------------------------------------------------------------------
#set up the python environment needed on lux13 for LDSC
# module purge
# module load miniconda/3.2021.10 ldsc/v1.0.1
# source activate ldsc

#set date 
DATE=`date +"%Y%m%d"`
#---------------------------------------------------------------------------------------------------
#set up our files: the general input path and the outputfolder
#and munged sumstats
#---------------------------------------------------------------------------------------------------
DTI_munged=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/DTI/LDSC/munged_enigmaEvo_DTI

#output gen cor
out_path_gencor=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/DTI/LDSC/gencor_phenoSPD
# details for LDSC
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB
ldscores_dir=${resource_dir}/LDscores

# #---------------------------------------------------------------------------------------------------
# #set up stuff for gen cor
# #----------------------------------------------------------------------------------------------------
# mkdir -p ${out_path_gencor}
# for mungedDTI_1 in ${DTI_munged}/*gz; do
	# filename_1=${mungedDTI_1##*/}
	# pheno_1=${filename_1##*europeanDTI_}
	# for mungedDTI_2 in ${DTI_munged}/*gz;do
		# filename_2=${mungedDTI_2##*/}
		# pheno_2=${filename_2##*europeanDTI_}
		# if [ ! -f ${out_path_gencor}/gencor_ukb43760_${pheno_1}_${pheno_2}.log ]; then
		# echo 'run LDSC to calculate genetic correlation between both traits' ${pheno_1} ${pheno_2}
		# ldsc.py \
		# --rg ${mungedDTI_1},${mungedDTI_2} \
		# --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
		# --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
		# --out ${out_path_gencor}/gencor_ukb43760_${pheno_1}_${pheno_2}
		# fi
	# done
# done


# ---------------------------------------------------------------------------------------------------
#extract gen cor resuls 
#----------------------------------------------------------------------------------------------------
#little wrapper to extract the summary values from the LDSC log after running gen cor
echo "trait rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se"> $out_path_gencor/gencor_DTI_phenoSPD.txt
for i in $out_path_gencor/*log; do
  res=$(cat $i | grep '^/d.*sumstat' | cut -d" " -f4-)
  pheno=$(basename "$i" | cut -d. -f1)
  echo $pheno $res >> $out_path_gencor/gencor_DTI_phenoSPD.txt
done 