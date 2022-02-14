#!/bin/bash
#$ -N AncestryRegression
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash

#-----Ancestry Regression-----
# Using a standard multivariable regression implemented in R with the lm() function, we regress 
# the GWAS summary statistics prior to ancestry regression (Beta_strat) with the ancestry 
# PCs (Beta_PCs). The residuals of this model are then ancestry corrected effect sizes (Beta_r) 
# and also have ancestry corrected standard errors and P-values. This is implemented in AncestryRegression_noGC.R.
# (taken from: https://bitbucket.org/jasonlouisstein/enigmaevolma6/src/master/1000Gphase3_PC_cor)

#-----
Rscript /data/clusterfs/lag/users/gokala/enigma-evol/AncestryRegression_noGC.R
