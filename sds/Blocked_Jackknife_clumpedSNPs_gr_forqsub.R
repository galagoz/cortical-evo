##/usr/local/R-3.2.3/bin/R
##Run Spearman's Correlation GWAS GenLang Data with SDS (~ 2000 years selection) using blocked jackknife
options(stringsAsFactors=FALSE)
library(GenomicRanges);

#load("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/ancreg/surfaceDK_parsorbitalis_globalCov_ancreg.Rdata")
#mergedGR=as.data.frame(mergedGR)
#write.csv(mergedGR,file="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/ancreg/surfaceDK_parsorbitalis_globalCov_ancreg.csv",row.names=F)
#csv_mergedGR=read.csv("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/ancreg/surfaceDK_parsorbitalis_globalCov_ancreg.csv")

##load SDS file
fSDS="/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/resources/SDS_UK10K_n3195_release_Sep_19_2016.tab.gz"
##directory of spearman's output
outputdir = "/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/results/replication/sds"
##Rdata files containing GWAS summary statistics
rdatafileloc = "/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface_ancreg/"
##read in gwas statistics file (compiled for all traits)
fGWASsumstats = "/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface_ancreg/sumstats_rdata_list.txt"

##Read in SDS file
SDS=read.table(fSDS, fill=TRUE, header=TRUE)

##Match the Rdata file locations of sumstats, text file sumstats, and clumped files
GWASsumstats=read.table(fGWASsumstats, header=FALSE)$V1;
##Parse to get trait name - UPDATE THIS PART ACCORDING TO YOUR PATH AND FILE NAMES
tmpname = sapply(GWASsumstats,function (x) {unlist(strsplit(as.character(x),"/",fixed=TRUE))[12]})
phenoname = sapply(tmpname,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[1]})
allfileloc = data.frame(rdatafile=GWASsumstats)

output=data.frame(global_corr_spearman=rep(NA, nrow(allfileloc)),
               BJK_ESTIM_AVE=rep(NA, nrow(allfileloc)),
               BJK_ESTIM_SE=rep(NA, nrow(allfileloc)),
               BJK_ESTIM_Z=rep(NA, nrow(allfileloc)),
               BJK_ESTIM_PVAL=rep(NA, nrow(allfileloc)))
#i=67
##Loop over each of the phenotypes
for (i in 1:nrow(allfileloc)) {
    pheno = phenoname[i]
    cat(' Working on:', pheno, '\n')
    cat('loading in',pheno,'pre-existing Rdata file...\n')
    #load(allfileloc$rdatafile[i])
    load(allfileloc$rdatafile[i])
    ##calculate Z score in GWAS files
    mergedGR$Z=NA
    mergedGR$Z=mergedGR$BETA/mergedGR$SE
    GWAS = as.data.frame(mergedGR) #removed mcols(mergedGR)
    ##Merge SDS with GWAS
    merged = merge(SDS, GWAS, by.x="ID", by.y="SNP") ##x=SDS, y=GWAS
    ##remove all NAs, keep only SNPs that have both measurements
    ind=which(!is.na(merged$Z)&!is.na(merged$SDS))
    merged=merged[ind,]
    ##for unaligned alleles, flip the SDS score
    unalignedind=which(toupper(merged$A1)==merged$AA & toupper(merged$A2)==merged$DA)

    merged$SDS[unalignedind] = (-1)* merged$SDS[unalignedind]
    tmpA1=merged$AA
    tmpA2=merged$DA
    merged$AA[unalignedind]=tmpA2[unalignedind]
    merged$DA[unalignedind]=tmpA1[unalignedind]
    ##Make sure alleles are aligned
    alignedind=which(toupper(merged$A1)==merged$DA & toupper(merged$A2)==merged$AA)
    merged=merged[alignedind,]
    ##Make sure to pick the trait increaing allele to change to tSDS
    negind=which(merged$Z < 0)
    tmpA1=merged$AA
    tmpA2=merged$DA
    merged$AA[negind]=tmpA2[negind]
    merged$DA[negind]=tmpA1[negind]
    merged$A1[negind]=tmpA1[negind]
    merged$A2[negind]=tmpA2[negind]
    merged$SDS[negind]= -1*merged$SDS[negind]
    merged$Z[negind]= -1* merged$Z[negind]

    ##Sort the merged file by genomic location
    newmergedGR = GRanges(merged$CHR,IRanges(merged$POS,merged$POS))
    mcols(newmergedGR) = merged[,c(1,4:17,25)] # for ENIGMA-replication ancreg sumstats, added 25th column (Z) here, otherwise "dat" didn't have newmergedGR$Z (19th column for non-ancreg sumstats)
    newmergedGR = sort(sortSeqlevels(newmergedGR))
    
    dat = cbind(newmergedGR$SDS, newmergedGR$Z)

    n_snps = dim(dat)[1]
    jackknife_num_blocks_flag = 100
    BJK_num_blocks=jackknife_num_blocks_flag
    BJK_block_size=floor(n_snps/BJK_num_blocks)
    BJK_num_blocks_p1=n_snps%%BJK_block_size

    ##tests: spearman corr > 0
    num_tests= 1
    
    ##Store leave-one-block-out estimates for different stats
    BJK_PSEUDO= matrix(0, nrow=BJK_num_blocks, ncol=num_tests)
    
    global_corr_spearman=cor(dat[,1],dat[,2],method="spearman")

    BJK_end=0
    BJK_start=1
    BJK_NatBrkPnt_start = 1
    BJK_NatBrkPnt_end   = 1

    for (BJK_iter in 1:BJK_num_blocks) {
    	##update
        BJK_start=BJK_end + 1

        if (BJK_iter <= BJK_num_blocks_p1){
           BJK_end = min(BJK_start + BJK_block_size,n_snps)
        } else {
           BJK_end = min(BJK_start + BJK_block_size - 1,n_snps)
        }

        if (BJK_start == 1) {
           working_indices = (BJK_end+1):n_snps
        } else if (BJK_end == n_snps) {
           working_indices = 1:(BJK_start-1)
        } else {
           working_indices = c(1:(BJK_start-1),(BJK_end+1):n_snps)
        }

      	## The Jackknife pseudovalues
        leave_block_out_corr_spearman=cor(dat[working_indices,1],dat[working_indices,2],method="spearman")
        BJK_PSEUDO[BJK_iter,1] = BJK_num_blocks * global_corr_spearman - (BJK_num_blocks-1) * leave_block_out_corr_spearman
   }

   ## Jackknife estimates of mean & standard-error => p-value

   BJK_ESTIM_AVE  = matrix(apply(BJK_PSEUDO,2,mean),nrow=1)
   BJK_ESTIM_SE   = matrix(apply(BJK_PSEUDO,2,sd)/sqrt(BJK_num_blocks),nrow=1)
   BJK_ESTIM_Z    = matrix(BJK_ESTIM_AVE/BJK_ESTIM_SE,nrow=1)
   BJK_ESTIM_PVAL = matrix(2*pnorm( -abs(BJK_ESTIM_Z) ),nrow=1)
   row.names(BJK_ESTIM_AVE) = "average"
   row.names(BJK_ESTIM_SE) = "error"
   row.names(BJK_ESTIM_Z) = "z-score"
   row.names(BJK_ESTIM_PVAL) = "p-value"
   BJK_ESTIM_AVE  =  signif(BJK_ESTIM_AVE,4)
   BJK_ESTIM_SE   =  signif(BJK_ESTIM_SE,4)
   BJK_ESTIM_Z    =  signif(BJK_ESTIM_Z,4)
   BJK_ESTIM_PVAL =  signif(BJK_ESTIM_PVAL,4)

   cat("Block-jackknife output\n")
   cat("Global Spearman corr: ", global_corr_spearman, "\n")

   output$global_corr_spearman[i]=global_corr_spearman
   output$BJK_ESTIM_AVE[i]=BJK_ESTIM_AVE
   output$BJK_ESTIM_SE[i]=BJK_ESTIM_SE
   output$BJK_ESTIM_Z[i]=BJK_ESTIM_Z
   output$BJK_ESTIM_PVAL[i]=BJK_ESTIM_PVAL

   rownames(output)[i]=pheno
   write.csv(output, file=paste0(outputdir, "/SDS_BJK_ancreg_replication_1kblocks.csv"))
}