#$ -N run_sumstats_Rdata2txt
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash

#-----Convert sumstats.Rdata to sumstats.txt-----

# Script to submit the conversion job to the grid. qsub this.
# Convert Rdata formatted summary statistics files to txt for downstream analysis.

Rscript sumstats_Rdata2txt.R
