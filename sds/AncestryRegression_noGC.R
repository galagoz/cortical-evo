##Run Correlation of effect sizes from GWAS of 1000G PC components with effect sizes from ENIGMA
##The effects sizes from GWAS of 1000G PC components were acquired from Katya and are from phase 3
##They were calculated by deriving PCs from 1000G (all populations) and correlating that with SNPs
##The goal here is to see if population stratification is driving the results
options(stringsAsFactors=FALSE);
library(GenomicRanges);

args = commandArgs(trailingOnly=TRUE)
##load 1000G PC effect sizes (basically an estimate population stratification for each SNP)
f1000G="/data/clusterfs/lag/users/gokala/genlang-evol/1kg_phase3_ns.allpop.unrel2261_eigenvec.P1to20_beta_se_pval.Rdata"
load(f1000G)
##directory of spearman's output
outputdir = "/data/clusterfs/lag/users/gokala/enigma-evol/ancreg/global/"
##Rdata files containing GWAS summary statistics
rdatafileloc = "/data/clusterfs/lag/users/gokala/enigma-evol/sumstatsRdata/global"
##read in gwas statistics file (compiled for all traits)
fGWASsumstats = "/data/clusterfs/lag/users/gokala/enigma-evol/sumstatsRdata/global/sumstats_rdata_list.txt"

##Match the Rdata file locations of sumstats, text file sumstats, and clumped files
GWASsumstats=read.table(fGWASsumstats, header=FALSE)$V1
##Parse to get trait name
tmpname = sapply(as.character(GWASsumstats),function (x) {unlist(strsplit(x,"/",fixed=TRUE))[10]})
phenoname = paste(sapply(tmpname,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[7]}),sapply(tmpname,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[8]}),sapply(tmpname,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[9]}),sep="_")
#phenoname = substr(tmpnamepheno,1,nchar(tmpnamepheno)-1)
allfileloc = data.frame(rdatafile=GWASsumstats)

#what are clumped files here?
output=data.frame(matrix(NA,ncol=41,nrow=nrow(allfileloc)))
colnames(output)[seq(1,40,2)] = paste0("clumpedcorestimate",seq(1,20))
colnames(output)[seq(2,40,2)] = paste0("clumpedp",seq(1,20))
colnames(output)[41] = "numclumpedsnps"
rownames(output) = phenoname

##Loop over each of the phenotypes
for (i in 1:nrow(allfileloc)) {
    pheno = phenoname[i];
    cat(' Working on:', pheno, '\n')
    if (!file.exists(allfileloc$rdatafile[i])) {
        cat('Reading Rdata file',pheno,'to save...\n');
        MA = read.table(allfileloc$txtfile[i],header=TRUE);
        ##Convert to genomic ranges and save as an Rdata file
        mergedGR = GRanges(MA$CHR,IRanges(MA$BP,MA$BP));
        mcols(mergedGR) = MA[,c(1:9)];
        save(mergedGR,file=allfileloc$rdatafile[i]);
    } else {
        cat('loading in',pheno,'pre-existing Rdata file...\n');
        load(allfileloc$rdatafile[i]);
    }
    #GWAS = as.data.frame(mcols(MAranges));
    ##Merge 1kgp with GWAS
    merged = merge(table, mergedGR, by="SNP") ##x=kG, y=GWAS
    ##remove all NAs, keep only SNPs that have both measurements
    ind=which(!is.na(merged$BETA) & !is.na(merged$P1beta)) # added P1beta NA filtering as well, to filter out SNPs that don't have PC loading values in 1kG file.
    merged=merged[ind,]

    ##for unaligned alleles, flip the Beta score
    unalignedind=which(toupper(merged$A1.x)==toupper(merged$A2.y) & toupper(merged$A2.x)==toupper(merged$A1.y));
    ##alignedind=which(toupper(merged$A1.x)==toupper(merged$A1.y) & toupper(merged$A2.x)==toupper(merged$A2.y));
    ##Flip the beta value to match the kG effect allele
    merged$BETA[unalignedind] = (-1)* merged$BETA[unalignedind];
    ##Remove the GWAS alleles to prevent confusion
    merged$A1.y = NULL;
    merged$A2.y = NULL;
    ##Run the regression to then remove the PC effects
    ##See https://www.biorxiv.org/content/biorxiv/early/2016/09/19/076133.full.pdf
    ##Correcting subtle stratification in summary association statistics
    modelfit = lm(merged$BETA ~ merged$P1beta + merged$P2beta + merged$P3beta + merged$P4beta + 
                      merged$P5beta + merged$P6beta + merged$P7beta + merged$P8beta + merged$P9beta + 
                      merged$P10beta + merged$P11beta + merged$P12beta + merged$P13beta + 
                      merged$P14beta + merged$P15beta + merged$P16beta + merged$P17beta + 
                      merged$P18beta + merged$P19beta + merged$P20beta); 
    ##Get the summary of the model
    smodelfit = summary(modelfit);
    #plot(merged$BETA,merged$P1beta)
    #abline(lm(merged$BETA ~ merged$P1beta,merged))
    ##These are the estimates of the ancestry corrected beta values
    ##Note that this does not include the intercept in the model (modelfit above)
    ##This makes it very slightly different than the resid(modelfit)
    betaPCs = smodelfit$coefficients[c(2:21),1];
    gammahat = merged[,seq(7,66,3)];
    ##Construct a yhat that is the linear combination of betaPCs and gammahat
    yhat = matrix(0,nrow=nrow(gammahat),ncol=1);
    for (iPC in 1:20) {
        yhat = yhat + betaPCs[iPC]*gammahat[,iPC];	
    }	
    merged$ancBETA = merged$BETA - yhat;
    merged$ancRESID = resid(modelfit)
    ##These are the estimates of the ancestry corrected SE values
    ##The formula here is from the paper above, slightly modified because it should be sqrt()
    ##The basics are that (when y and yhat are independent), var(y-yhat) = var(y) + var(yhat)
    varyhat = 0;
    segammahat = merged[,seq(8,66,3)];
    for (iPC in 1:20) {
        varyhat = varyhat + (betaPCs[iPC])^2*(segammahat[,iPC])^2;	
    }
    merged$ancSE = sqrt( (merged$SE)^2 + varyhat);
    ##Finally turn this into a P-value
    ##Paper says that ancBETA^2 / ancSE^2 is distributed as chisq with 1 d.f.
    merged$ancP = pchisq((merged$ancBETA)^2/(merged$ancSE)^2,1,lower.tail=FALSE);
    ##Convert back to Genomic Ranges and save the data
    mergedGR = GRanges(merged$CHR.x,IRanges(merged$POS,merged$POS));    
    mcols(mergedGR) = merged[,c(1,5,6,67:76)];
    mergedGR = sort(sortSeqlevels(mergedGR));
    save(mergedGR,file=paste0(outputdir,pheno,"_ancreg.Rdata"))
    ##Merge with the clumped list to get only LD independent SNPs
    ##clumpedlist=read.table(allfileloc[i,3],fill=TRUE,header=TRUE);
    ##clumpedmerged=merge(merged, clumpedlist, by="SNP");
    ##output$numclumpedsnps[i] = nrow(clumpedmerged);
    ##Loop over all components to determine residual correlation with SNP-PC effect sizes
    ##for  (k in 1:20) {  
    ##	 corout.clumped=cor.test(clumpedmerged[,seq(8,66,3)[k]], clumpedmerged$ancBETA, method="pearson");
    ##	 output[i,seq(1,40,2)[k]]=corout.clumped$estimate
    ##	 output[i,seq(2,40,2)[k]]=corout.clumped$p.value
    ##}
    ##write.csv(output, paste0(outputdir, "/corvalues_residualCLUMP1.csv"))
}
