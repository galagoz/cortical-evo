#!/usr/bin/env Rscript
#
# This script will parallelize Overlap Analysis 
# jobs over annotations.
# The R script works only with R version >4.0.3.
# Thus, run this script on gridportal1 instead of
# gridmaster.
#
# Gokberk Alagoz
# Created on: 02.07.2021
#--------------------------------------------------#
# VARIABLES

bedfilesDir=$1
sumstatsList=$2
clumpedSumstatsDir=$3
outDir=$4

#--------------------------------------------------#

mkdir "${outDir}/scripts"
mkdir "${outDir}/results"
mkdir "${outDir}/logs"

#--------------------------------------------------#
# FUNCTIONS

for i in $bedfilesDir/*.bed; do
	
	annot=$(basename "$i")
	shellFile="${outDir}/scripts/${annot}.sh"
	logFile="${outDir}/logs/${annot}.out"

	echo "Submitting jobs for " $annot

	echo '#$ -N eqtl_overlap
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash

module load R/R-4.0.3
Rscript /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/scripts/enigmaevol/eqtl/overlap_analysis_eQTL.R '$i' '$sumstatsList' '$clumpedSumstatsDir' '${outDir}/results'' > $shellFile

	chmod a+x ${shellFile}
	echo "Submitting $annot to the Grid"
	qsub -o ${logFile} -j y ${shellFile}

done < $sumstatsList

echo "All submitted!"

#--------------------------------------------------#
