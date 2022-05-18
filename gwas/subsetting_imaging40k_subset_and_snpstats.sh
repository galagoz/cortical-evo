#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -m eas
#$ -N subset_snpstats_DTI
#$ -q single15.q
#$ -e /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/log
#$ -o /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/log
#$ -M Barbara.Molz@mpi.nl

#script written by Dick Schijven as part of the UKBB Imaging preprocessing

usage()
{
   echo -e "\n-------------------------------------------------------------------------------------------------------------------------------------------------------------"
   echo "SUBSETTING AND SNP STATISTICS CALCULATION FOR THE IMAGING 40K DATASET"
   echo " "
   echo "This script extracts a custom list of subjects from the full 40k imaging genetic dataset. In addition, a .bgen index file and .sample file are generated"
   echo "and SNP statistics are calculated."
   echo " "
   echo "Usage: bash imaging40k_subset_and_snpstats.sh -s [configfile] -c [chr]"
   echo -e "\t-c [chr] Chromosome number. Possible values are 1-22 and X"
   echo -e "\t-s [configfile] Configuration file for the run. Use the general template."
   echo -e "\t-h Show help."
   echo " "
   echo "When running the script on the cluster (through gridmaster), you might want to provide arguments specific to qsub, such as the location where qsub stores a"
   echo "standard log file or the email address to which a message should be sent when the script finished, in addition to the arguments specific to the script."
   echo "If this is the case, provide the qsub arguments before the script name, and the arguments specific to the script after the script name:"
   echo "qsub -M Your.Name@mpi.nl -e /dir/where/err/is/saved -o dir/where/log/is/saved imaging40k_subset_and_snpstats.sh -s [configfile] -c [chr]"
   echo -e "-------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
   exit 1 # Exit script after printing usage
}

while getopts "hc:s:"  opt
do
   case "$opt" in
      c ) chr="$OPTARG" ;;
      s ) configfile="$OPTARG" ;;
	  h ) usage
	  	  exit 1 ;;
      ? ) usage ;; # Print usage in case parameter is non-existent
   esac
done

# If one of the required parameters, chr of configfile, is empty, show an error message
if [ -z ${chr} ] || [ -z ${configfile} ]; then
	echo -e "\nERROR: Required input arguments are missing. Specify a chromosome number and the path of the configuration file when running the script.";
	usage
fi

# Set directory where qctool version 2.0.6 is located
qctoolDir=/data/workspaces/lag/shared_spaces/Resource_DB/qctool_v2.0.6/

# Set the directory where bgenix software is located
bgenixDir=/data/workspaces/lag/shared_spaces/Resource_DB/bgen/gavinband-bgen-af23766eee3d/build/apps/

# Read variables from the configuration file
source ${configfile}

# Print a standard header for the log file
printf "IMAGING 40K DATASET SUBSETTING FOR CHROMOSOME "${chr}"\n\n"
printf "Run name:\t\t"${RunName}"\n"
printf "Configuration file:\t"${configfile}"\n"

printf "Input .gen file:\t"${ImagingSubsetData}"_chr"${chr}".bgen\n"
printf "Input .sample file:\t"${ImagingSubsetData}"_chr"${chr}".sample\n\n\n"


##GENERATE BGEN FILES

if [ -f ${OutDir}/${RunName}_chr${chr}.bgen ]; then

	printf "BGEN file for chr "${chr}" already exists.\n\n"

elif [ ${chr} == "X" ]; then

	printf "Subset BGEN file for chr X, including samples in "${SubjectList}"\n\n"
	${qctoolDir}/qctool \
	-g ${ImagingSubsetData}_chr${chr}.bgen \
	-s ${ImagingSubsetData}_chr${chr}.sample \
	-incl-samples ${SubjectList} \
	-ofiletype "bgen_v1.2" \
	-og ${OutDir}/${RunName}_chr${chr}.bgen \
	-log ${OutDir}/${RunName}_chr${chr}.bgen.log

elif [ ${chr} -ge 1 ] && [ ${chr} -le 22 ]; then
 	
	printf "Subset BGEN file for chr "${chr}", including samples in "${SubjectList}"\n\n"
	${qctoolDir}/qctool \
	-g ${ImagingSubsetData}_chr${chr}.bgen \
	-s ${ImagingSubsetData}_chr${chr}.sample \
	-incl-samples ${SubjectList} \
	-ofiletype "bgen_v1.2" \
	-og ${OutDir}/${RunName}_chr${chr}.bgen \
	-log ${OutDir}/${RunName}_chr${chr}.bgen.log

else

	printf "\nERROR: Input chromosome number incorrect. Possible values for chromosome are 1-22 or X.\n"
	
	exit 1

fi


##GENERATE BGEN INDEX FILE

if [ -f ${OutDir}/${RunName}_chr${chr}.bgen.bgi ]; then

	printf "BGEN index file for chr "${chr}" already exists.\n\n"

elif [ ${chr} -ge 1 ] && [ ${chr} -le 22 ] || [ ${chr} == "X" ]; then

	${bgenixDir}/bgenix -g ${OutDir}/${RunName}_chr${chr}.bgen -index

else

	printf "\nERROR: Input chromosome number incorrect. Possible values for chromosome are 1-22 or X.\n"
	
	exit 1

fi


##GENERATE SAMPLE FILES

if [ -f ${OutDir}/${RunName}_chr${chr}.sample ]; then

	printf "Sample file for chr "${chr}" already exists.\n\n"

elif [ ${chr} == "X" ]; then
 	
 	printf "Writing sample file for "${OutDir}"/"${RunName}"_"${chr}".bgen\n\n"
	
	${qctoolDir}/qctool \
	-g ${OutDir}/${RunName}_chr${chr}.bgen \
	-os ${OutDir}/${RunName}_chr${chr}.sample \
	-log ${OutDir}/${RunName}_chr${chr}.sample.log

	mv ${OutDir}/${RunName}_chr${chr}.sample ${OutDir}/${RunName}_chr${chr}_.sample
	awk '(NR==FNR){arr[$1];next}($1 in arr)' ${OutDir}/${RunName}_chr${chr}_.sample ${ImagingSubsetData}_chr${chr}.sample > ${OutDir}/${RunName}_chr${chr}.sample
	paste ${OutDir}/${RunName}_chr${chr}_.sample ${OutDir}/${RunName}_chr${chr}.sample > ${OutDir}/${RunName}_compareX.sample
	awk '$1 != $2 {print $0}' ${OutDir}/${RunName}_compareX.sample > ${OutDir}/${RunName}_compareX.nonmatching.sample

	if [ $( wc -l ${OutDir}/${RunName}_compareX.nonmatching.sample | awk '{print $1}' ) -gt 1 ]; then

		printf "\nERROR: Subject list mismatch when constructing the .sample file for chromosome X.\n"
		exit 1

	fi

	rm ${OutDir}/${RunName}_chrX_.sample; rm ${OutDir}/${RunName}_compareX.sample; rm ${OutDir}/${RunName}_compareX.nonmatching.sample

elif [ ${chr} -ge 1 ] && [ ${chr} -le 22 ]; then

	printf "Writing sample file for "${OutDir}"/"${RunName}"_"${chr}".bgen\n\n"
	
	${qctoolDir}/qctool \
	-g ${OutDir}/${RunName}_chr${chr}.bgen \
	-os ${OutDir}/${RunName}_chr${chr}.sample \
	-log ${OutDir}/${RunName}_chr${chr}.sample.log

fi


##GENERATE SNP-STATS FILES

if [ -f ${OutDir}/${RunName}_chr${chr}.snpstats.txt ]; then

	printf "SNP-stats file for chr "${chr}" already exists.\n\n"

elif [ ${chr} == "X" ]; then

 	printf "Writing snp-stats file for "${OutDir}"/"${RunName}"_"${chr}".bgen\n\n"

 	${qctoolDir}/qctool \
	-g ${OutDir}/${RunName}_chr${chr}.bgen \
	-s ${OutDir}/${RunName}_chr${chr}.sample \
	-infer-ploidy-from sex \
	-snp-stats \
	-osnp ${OutDir}/${RunName}_chr${chr}.snpstats.txt \
	-log ${OutDir}/${RunName}_chr${chr}.snpstats.log

elif [ ${chr} -ge 1 ] && [ ${chr} -le 22 ]; then
 	
 	printf "Writing snp-stats file for "${OutDir}"/"${RunName}"_"${chr}".bgen\n\n"
	
	${qctoolDir}/qctool \
	-g ${OutDir}/${RunName}_chr${chr}.bgen \
	-snp-stats \
	-osnp ${OutDir}/${RunName}_chr${chr}.snpstats.txt \
	-log ${OutDir}/${RunName}_chr${chr}.snpstats.log


else

	printf "\nERROR: Input chromosome number incorrect. Possible values for chromosome are 1-22 or X.\n"
	
	exit 1

fi
