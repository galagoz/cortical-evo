##/usr/local/R-3.2.3/bin/R
##Run Correlation of effect sizes from GWAS of 1000G PC components with effect sizes from ENIGMA
##The effects sizes from GWAS of 1000G PC components were acquired from Katya and are from phase 3
##They were calculated by deriving PCs from 1000G (all populations) and correlating that with SNPs
##The goal here is to see if population stratification is driving the results

##This script does the work for run_PC_cor_BJK_noGC.R
 
options(stringsAsFactors=FALSE)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)
frdatafile = args[1];
phenoname = args[2];
outputdir = args[3];

print(frdatafile)
print(phenoname)
print(outputdir)

##frdatafile = "Rdatafiles_noGC/Mean_bankssts_surfavg.Rdata";
##phenoname = "Mean_bankssts_surfavg";
##outputdir = "/ifs/loni/faculty/dhibar/ENIGMA3/MA6/evolution/1000Gphase3_PC_cor/beta_strat_noGC/output";

##load 1000G PC effect sizes (basically an estimate population stratification for each SNP)
f1000G="/data/clusterfs/lag/users/gokala/genlang-evol/1kg_phase3_ns.allpop.unrel2261_eigenvec.P1to20_beta_se_pval.Rdata"

##Block jackknife function
cor.test.jack = function(x, y, blocks=1000) {

  keep = is.finite(x) & is.finite(y)
  x = x[keep]
  y = y[keep]
  mat = cbind(x, y, floor(seq(1, blocks+1-1e-3, len=length(x))))
  part_ests = sapply(1:blocks, function(i) cor(mat[mat[,3] != i, 1], mat[mat[,3] != i, 2]))
  jack_est = mean(part_ests)
  jack_se = sqrt((blocks-1)/blocks * sum((jack_est - part_ests)^2))
  jack_z = jack_est/jack_se;
  jack_p = 2*pnorm( -abs(jack_z));
  list(estimate=jack_est, se=jack_se, p=jack_p, estimate_normal=cor(x, y), se_normal=cor.test.plus(x, y)$se)

}

cor.test.plus <- function(x, y, ...) {

  # like cor.test, but also returns se of correlation
  corr = cor.test(x, y, ...)
  corr$se = unname(sqrt((1 - corr$estimate^2)/corr$parameter))
  corr

}

##set output variable
output=data.frame(matrix(NA,ncol=60,nrow=1));
colnames(output)[seq(1,60,3)] = paste0("BJK_cor",seq(1,20));
colnames(output)[seq(2,60,3)] = paste0("BJK_SE",seq(1,20));
colnames(output)[seq(3,60,3)] = paste0("BJK_P",seq(1,20))
rownames(output) = phenoname;

##Read in 1000G file
load(f1000G);
##Read in the summary statistics
load(frdatafile);

##calculate Z score in GWAS files
mergedGR$Z=NA;
mergedGR$Z=mergedGR$BETA/mergedGR$SE;
GWAS = mergedGR ### THIS LINE WAS ADDED
#GWAS = as.data.frame(mcols(tmp_ss_table)); ### THIS LINE IS REMOVED 

##Merge 1kgp with GWAS
merged = merge(table, GWAS, by="SNP") ##x=kG, y=GWAS
##remove all NAs, keep only SNPs that have both measurements
ind=which(!is.na(merged$Z))
merged=merged[ind,]
##for unaligned alleles, flip the Z score
unalignedind=which(toupper(merged$A1.x)==toupper(merged$A2.y) & toupper(merged$A2.x)==toupper(merged$A1.y));
##alignedind=which(toupper(merged$A1.x)==toupper(merged$A1.y) & toupper(merged$A2.x)==toupper(merged$A2.y));

merged$Z[unalignedind] = (-1)* merged$Z[unalignedind];
merged$BETA[unalignedind] = (-1)* merged$BETA[unalignedind];

##Need to resort the merged file
mergedGR = GRanges(merged$CHR.x,IRanges(merged$POS,merged$POS))
mcols(mergedGR) = merged[,c(1,3,5:74)]
merged = as.data.frame(sort(sortSeqlevels(mergedGR)))

##Loop over all components
for  (k in 1:20) {

     coroutput=cor.test.jack(merged[,seq(10,67,3)[k]], merged$BETA)
     output[1,seq(1,60,3)[k]] = coroutput$estimate
     output[1,seq(2,60,3)[k]] = coroutput$se
     output[1,seq(3,60,3)[k]] = coroutput$p

     print(k)
     print(coroutput)

}
  
write.csv(output, paste0(outputdir, "/corvalues_",phenoname,"_BJK.csv"))
