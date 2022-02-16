#$ -N sds_ancreg
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash

#-----Run SDS-----
# Master script to run SDS on the cluster. qsub this.

#-----
Rscript Blocked_Jackknife_clumpedSNPs_gr_forqsub.R
#-----
