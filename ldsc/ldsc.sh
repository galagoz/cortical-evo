#!/bin/sh
#-----Run ldsc-----

# Estimating heritability, genetic correlation and the LD Score Regression Intercept
# ldsc function is from github.com/bulik/ldsc
# For interpreting results, check the tutorial at github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation

#-----Variables-----
# $ss1 - munged summary statistics file
# $ss2 - ancestry regressed + munged summary statistics file
# $iv - LD Score files to use as the independent variable in the LD Score regression
# $rw - LD Scores to use for the regression weights
# For this analysis, same file is used as iv and rw -> eur_w_ld_chr

inDir="/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface/withGlob/munged/"
outDir="/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface/withGlob/munged/ldsc/"
iv_rw="/data/workspaces/lag/shared_spaces/Resource_DB/LDscores/eur_w_ld_chr/"
mungedSumstatsList="/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface/withGlob/munged/munged_sumstats_list.txt"

#-----
mkdir "${inDir}scripts"
mkdir "${inDir}ldsc"

echo $inDir
echo $outDir
echo $iv_rw
echo "Starting to compute LDSCs"

#module load python/3.7.2 \
#module load ldsc/v1.0.1 \

while read line; do
   echo $line
   tmp_base_name=$(basename "$line" .txt)
   echo $tmp_base_name
   pheno_name="$(cut -d'_' -f5- <<<"$tmp_base_name")"
   echo $pheno_name
   output="${outDir}${pheno_name}"
   tmp_run_file="${inDir}scripts/${pheno_name}.sh"
   echo '#!/bin/sh
#$ -N ldsc
#$ -cwd
#$ -q single.q
#$ -S /bin/bash

/home/gokala/programs/ldsc/ldsc.py --h2 '$line' --ref-ld-chr '$iv_rw' --w-ld-chr '$iv_rw' --out '$output'' > $tmp_run_file
   chmod a+x $tmp_run_file
   echo "Created the script for cluster ->  submitting ${pheno_name} to the Grid"
   qsub -wd "${inDir}scripts" $tmp_run_file
done < $mungedSumstatsList

echo "Done!"

#-----
