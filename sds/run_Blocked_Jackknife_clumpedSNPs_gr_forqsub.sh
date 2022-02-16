#!/bin/sh
#$ -N sds_BJK
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash

#-----Run SDS-----

#

#-----
Rscript Blocked_Jackknife_clumpedSNPs_gr_forqsub.R
#-----
