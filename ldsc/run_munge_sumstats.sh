#!/bin/sh
# Gokberk Alagoz, 2020

#-----Munge Sumstats-----

# Reformat GWAS summary stats before computing LDSC intercept
# munge_sumstats.py is from github.com/bulik/ldsc
# Based on the tutorial at github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
# WARNING: Run it on single.q, other queues won't work because of Python and Numpy
# version incompatibility.

#-----Variables-----
# $inDir - path to summary statistics files
# $outDir - path to your output folder
# $sumstatsList - a text file with txt formatted sumstats files with full paths

inDir="/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface_ancreg/"
outDir="/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface_ancreg/munged/"
sumstatsList="/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface_ancreg/sumstats_txt_list.txt"

#-----

echo $inDir
echo $outDir
echo "Munging sumstats"
mkdir "${inDir}scripts"
mkdir "${inDir}munged"

while read line; do
   echo $line
   LINE=$line
   tmp_base_name=$(basename "$line" .txt)
   echo $tmp_base_name
   pheno_name="$(cut -d'_' -f1- <<<"$tmp_base_name")"
   echo $pheno_name
   tmp_run_file="${inDir}scripts/${pheno_name}.sh"
   output="${outDir}${tmp_base_name%.txt}_munged.txt"
   echo $line
   echo '#!/bin/bash
#$ -N munge_sumstats
#$ -cwd
#$ -q single.q
#$ -S /bin/bash

python /home/gokala/programs/ldsc/munge_sumstats.py --sumstats '$LINE' --out '$output' --merge-alleles /data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/resources/w_hm3.snplist' > $tmp_run_file
   chmod a+x $tmp_run_file
   echo "Created the script for cluster ->  submitting ${pheno_name} to the Grid"
   qsub -wd "${inDir}scripts" $tmp_run_file
done < $sumstatsList

echo "Finished!"

#-----
