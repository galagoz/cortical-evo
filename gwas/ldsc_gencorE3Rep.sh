
#----------------------------------------------------------------
#This scrip runs LDSC genetic correlations between derived summary statistics and Enigma 3 summary statistics 
#Info about LDSC: https://github.com/bulik/ldsc
#Info about summary stats from E3: https://www.science.org/doi/10.1126/science.aay6690
# by Barbara Molz, Nov 2021


#!/bin/sh
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
ukb_munged=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/replication/LDSCUpdate/mungedReplication
e3_sumstats=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/e3_sumstats

#output gen cor
out_path_gencor=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/replication/LDSCUpdate/gencor_E3vsRep
# details for LDSC
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB
ldscores_dir=${resource_dir}/LDscores

#---------------------------------------------------------------------------------------------------
#set up stuff for gen cor
#----------------------------------------------------------------------------------------------------
mkdir -p ${out_path_gencor}

for sumstats in ${ukb_munged}/*gz; do
	if [[ $sumstats == *"Full"* ]]; then
		filename=${sumstats##*/}
		pheno=`echo ${filename} | awk -F_ '{print $4}'`
		e3_file=$e3_sumstats/enigma3_global_${pheno}.sumstats.gz
		if [ ! -f ${outpat}.log ]; then
			echo 'run LDSC to calculate genetic correlation between E3 and enigma Evo for global' ${pheno}
			ldsc.py \
			--rg ${e3_file},${sumstats} \
			--ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
			--w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
			--out $out_path_gencor/gencorE3_Rep_global_${pheno}
		fi
	else
		filename=${sumstats##*/}
		pheno=`echo ${filename} | awk -F_ '{print $4}'`
		roi=`echo ${filename} | awk -F_ '{print $5}'`
		covar=`echo ${filename} | awk -F_ '{print $6}'`
		e3_file=$e3_sumstats/enigma3_${roi}_${pheno}_${covar}.sumstats.gz
		if [ ! -f ${outpat}.log ]; then
			echo 'run LDSC to calculate genetic correlation between E3 and enigma Evo for regional' ${pheno} ${roi} ${covar}
			ldsc.py \
			--rg ${e3_file},${sumstats} \
			--ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
			--w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
			--out $out_path_gencor/gencorE3_Rep_${pheno}_${roi}_${covar}
		fi
	fi
done


#------------------------------------
# create summary table
#------------------------------------
#little wrapper to extract the summary values from the LDSC log after running gen cor
echo "trait rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se"> $out_path_gencor/summary_gencor_E3_${sample}.txt
for i in $out_path_gencor/*log; do
  res=$(cat $i | grep '^/d.*sumstat' | cut -d" " -f4-)
  pheno=$(basename "$i" | cut -d. -f1)
  echo $pheno $res >> $out_path_gencor/summary_gencor_E3_${sample}.txt
done 

