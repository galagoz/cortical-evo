#$ -N convert_txt2Rdata
#$ -cwd
#$ -q single.q
#$ -S /bin/bash

#-----Convert sumstats.txt to sumstats.Rdata-----

# Convert txt formatted summary statistics files to Rdata for downstream analysis.

#-----



echo "Starting to convert sumstats from txt to Rdata"
/usr/local/bin/Rscript sumstats_txt2Rdata.R
echo "Done!"




#-----
