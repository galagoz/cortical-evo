#-----Run Partitioned Heritability-----

# Estimating heritability, genetic correlation and the LD Score Regression Intercept
# ldsc function is from github.com/bulik/ldsc
# For interpreting results, check the tutorial at github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
# This script takes 2 parameters: A csv file with "Annotation" and "Baseline" columns and a text file list of munged sumstats.

#-----Variables-----
# $inDir - ancestry regressed + munged summary statistics directory

annotFile=$1
mungedSumstats=$2
firstLine=$(head -n 1 $mungedSumstats)
inDir=$(dirname $firstLine)
outDir="${inDir}/results/"

#-----
mkdir "${inDir}/scripts"
mkdir "${inDir}/results"
mkdir "${inDir}/shelloutput/"

echo $inDir
echo $outDir
echo "Starting to compute partitioned heritability"

while read a; do
   echo $a
   tmp_file=$(basename "$a" .txt)
   echo $tmp_file
   tmp_pheno="$(cut -d'_' -f4- <<<"$tmp_file")"
   echo $tmp_pheno
   #tmp_output="${outDir}${annot}"
   
   sed 1d $annotFile | while IFS="	" read -r col1 col2; do

	   annot=$col1
	   baseline=$col2
	   shellFile="${inDir}/scripts/${annot}_${tmp_pheno}.sh"
	   logFile="${inDir}/shelloutput/${annot}_${tmp_pheno}.out"
	   mkdir "${outDir}${annot}"
	   tmp_output="${outDir}${annot}/${tmp_file}"

	   echo '************************************************'
	   echo $tmp_pheno
	   echo $annot
	   echo $baseline
	   echo '************************************************'
	   
	   if [ $baseline == "baseline" ]
	   then
	   
	   echo '#$ -N partherit
#$ -cwd
#$ -q single.q
#$ -S /bin/bash

/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/scripts/cortical-evo/partherit/analysis/partherit_baseline_replication.sh '$a' '$annot' '$tmp_output'' > $shellFile
   	
	   chmod a+x ${shellFile}
	   echo "Created the script for cluster ->  submitting ${tmp_pheno} to the Grid"
	   qsub -o ${logFile} -j y ${shellFile}
	   
   	   else
	   
	   echo '#$ -N partherit
#$ -cwd
#$ -q single.q
#$ -S /bin/bash

/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/scripts/cortical-evo/partherit/analysis/partherit_extraAnnot_replication.sh '$a' '$annot' '$tmp_output' '$baseline'' > $shellFile
	   
	   chmod a+x ${shellFile}
	   echo "Created the script for cluster ->  submitting ${tmp_pheno} to the Grid"
	   qsub -o ${logFile} -j y ${shellFile}

	   fi

   done;

done < $mungedSumstats

echo "Done!"

#-----
