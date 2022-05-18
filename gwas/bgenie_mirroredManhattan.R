#----------------------------------------------------------------------------------------------------------
# This script will use new R package to plot mirrored manhattan plots to visualize Repilcation vs E3 or Le/Re'
# by Barbara Molz January 2022
#----------------------------------------------------------------------------------------------------------

# 
if("readxl" %in% rownames(installed.packages()) == FALSE) {install.packages("readxl")}
#check if packages exist and install in case
library("data.table")
library(readxl)
#this is the library needed for manhattan mirror plots. Details see: https://github.com/anastasia-lucas/hudson/
library(hudson)

#set paths to GWAS sumstats
setwd(paste0('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final'))
#gwas_left = fread(paste0('./european/bgenieOutput/global/sumstats_ukb43760_global_surface_le_european_allChr.gz'), header = T, select= c(1,10,11,7))
#gwas_right= fread(paste0('./european/bgenieOutput/global/sumstats_ukb43760_global_surface_re_european_allChr.gz'), header = T, select= c(1,10,11,7))
gwas_replication=fread('./replication/bgenieOutput/global/sumstats_ukb43760_global_surface_replication_allChr.gz', header = T, select= c(1,10,11,7))
gwas_e3=fread('ENIGMA3_mixed_se_wo_Mean_Full_SurfArea_20190429_noGC.txt.gz', header = T, select= c(1,10,11,7))

#change columnames to  ones recognized by the program
colnames(gwas_replication)[4] <- "pvalue"
colnames(gwas_replication)[3] <- "POS"

colnames(gwas_e3)[4] <- "pvalue"
colnames(gwas_e3)[3] <- "POS"

#colnames(gwas_left)[4] <- "pvalue"
#colnames(gwas_left)[3] <- "POS"

#colnames(gwas_right)[4] <- "pvalue"
#colnames(gwas_right)[3] <- "POS"

#only get the columns we need -
gwas_e3 <- gwas_e3[, c(1, 3, 4, 2)]
gwas_replication <- gwas_replication[, c(1, 3, 4, 2)]
png('./plots/miamiPlot_replicaton_E3.png', width = 1000, type = "cairo-png", res=300)
gmirror(top=gwas_replication, bottom=gwas_e3, chroms = c(1:22),tline=8.3e-10, bline=8.3e-10, 
        toptitle="cortical surface: replication", bottomtitle = "cortical surface: Enigma 3",file="miamiPlot_replicaton_E3",hgt=10)
dev.off()

#png(paste0('./plots/manhattan_twin_left_right.png'),width = 1000, type = "cairo")
#gmirror(top=gwas_left, bottom=gwas_right,chroms = c(1:22),tline=8.3e-10, bline=8.3e-10, 
 #       toptitle="cortical surface: left", bottomtitle = "cortical surface: right")
#dev.off()