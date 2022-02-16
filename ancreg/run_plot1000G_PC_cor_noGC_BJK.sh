#!/bin/bash
#$ -N wo_ancreg_plot
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash

#-----Plot Correlation Test Results (1000G PC loadings vs. sum.stats. wo/ancestry regression)-----
# 

#-----
Rscript /data/clusterfs/lag/users/gokala/enigma-evol/plot1000G_PC_cor_noGC_BJK.R
