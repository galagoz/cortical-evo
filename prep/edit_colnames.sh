#!/bin/bash

# Gokberk Alagoz

# This script will...
# After converting ancestry regressed Rdata sumstats to txt format,
# A1.x, A2.x, CHR.y colnames must be corrected before munging.
# Also, ancBETA column should be moved to 5th field and overwrite BETA column,
# so that ancestry regressed BETA values are used for LDSC partherit analysis.

dir="/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface_ancreg"

for file in $dir/*ancreg.txt; do

        awk '{print $1, $2, $3, $4, $13, $6, $7, $8, $9, $10, $11, $12}' $file > ${file%.txt}_formatted.txt
	sed -Ei '1s/A1.x/A1/;1s/A2.x/A2/;1s/CHR.y/CHR/;1s/ancBETA/BETA/' ${file%.txt}_formatted.txt;

done
