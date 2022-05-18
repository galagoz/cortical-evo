#----------------------------------------------------------------------------------------------------------
# script for the calculation of the number of independent phenotypes in GenLang
# by Else Eising, 11 June 2019
# adabted by Barbara Molz, May 2022 
# The script uses PhenoSpD
# see https://github.com/MRCIEU/PhenoSpD and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6109640/#sup6
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
# Set directories
#----------------------------------------------------------------------------------------------------------

PhenospdDir="/home/barmol/PhenoSpD/"
LDSCdir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/LDSCUpdate/gencor_phenoSPD/genCor_surfaceHemi.txt"
PhenospdOuptupDir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/LDSCUpdate/PhenoSpD/"

#----------------------------------------------------------------------------------------------------------
# Install and update PhenoSpD
#----------------------------------------------------------------------------------------------------------

# mkdir $PhenospdDir
# cd /home/barmol/
# git clone https://github.com/MRCIEU/PhenoSpD.git
# cd $PhenospdDir
# git pull


#----------------------------------------------------------------------------------------------------------
# Run PhenoSpD on LDSC correlation matrix
#----------------------------------------------------------------------------------------------------------


mkdir -p ${PhenospdOuptupDir}
cd $PhenospdDir
##using the created UKBiobank Genetic correlations correlation
Rscript ./script/phenospd.r --phenocorr $LDSCdir --out $PhenospdOuptupDir/GenCor_SurfacebasedHemi_phenoSPD



