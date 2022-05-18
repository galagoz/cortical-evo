#!/bin/sh
#$ -cwd
#$ -q single.q
#$ -e /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/log
#$ -o /data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/log
#$ -S /bin/bash
#$ -m eas
#$ -M Barbara.Molz@mpi.nl

chr=$1
pheno=$2



#run variant_qc
config_file=vqc_imaging40k_variant_qc_config_${pheno}_ukb43760.R

printf "RUN VARIANT-LEVEL QC with ${config_file} FOR CHROMOSOME "${chr}"\n\n"
Rscript vqc_imaging40k_variant_qc.R ${config_file} ${chr}

printf "\n\n"
