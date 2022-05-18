#----------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------
# This script uses the ouptut created by begenie_merge.py to create manhatten and qq plots for all DK hemiphsphere (total surface area/thickness)
# by Barbara Molz January 2022
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#load packages 
if("qqman" %in% rownames(installed.packages()) == FALSE) {install.packages("qqman")}
#check if packages exist and install in case
if("readxl" %in% rownames(installed.packages()) == FALSE) {install.packages("readxl")}
#check if packages exist and install in case
library(qqman)
library("data.table")
library(readxl)

sample_set='european'
#set up wd and create empty list that can hold the lambda details
setwd(paste0('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/',sample_set,'/le_re/sumstats/'))
lambdalist = list()
#start counter for later indexing
count=1
#import the roi list 
roi_list = read_excel(paste0('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/QC_5MAD/ROI_table.xlsx'),col_names = 'ROI')
#-------------------------------------------------------------------------------------
# loop through phenotypes, load one DK ROI at a time
for (pheno in c('surface', 'thickness'))
{
  for (hemi in c('le', 're'))
  {
    for (roi in roi_list$ROI)
    {
      for (covar in c('withGlob', 'withoutGlob'))
      {
        # only load stuff we need which is Chr, pos, rsid and one phenotype --> NOTE they have been converted back to pvalues during output merge
        gwas = fread(paste0('sumstats_ukb43760_regionalDK_',pheno,'_',hemi,'_',roi,'_', covar,'_',sample_set,'_allChr.gz'), header = T, select= c(1,7,10,11))
        # get the phenotype name to create a plot title
        cur_pheno = paste0('ukb43769_regionalDK_',pheno,'_',hemi,'_',roi,'_',covar,'_', sample_set)
        #save the png
        png(paste0('plots/bgenie_manhattan_',cur_pheno,'_.png'),width = 1000, type = "cairo")
        #create the manhattan with certain settings
        pmin=ceiling(abs(log10(min(gwas$P))))
        manhattan(gwas,chr="CHR",bp ="BP",p = 'P',snp= 'SNP', genomewideline =-log10(8.3e-10) , suggestiveline=-log10(5e-08), ylim=c(0,pmin), chrlabs = as.character(1:22), main = cur_pheno)
        dev.off()
        chisq <- qchisq(1-gwas$P,1)
        mchisq<-summary(chisq)[c("Mean","Median","Max.")]
        ## calculate lambda gc
        lambda <- round(median(chisq,na.rm=TRUE)/qchisq(0.5,1),digits=5)
        # save lambda values into file
        t<-cbind(name=cur_pheno,lambda_all=lambda,mean_chisq_all=mchisq[1],median_chisq_all=mchisq[2],max_chisq_all=mchisq[3])
        t<-as.data.frame(t)
        lambdalist[[count]] <-t
        png(paste0('plots/bgenie_qqplot_',cur_pheno,'_.png'),width = 1000, type = "cairo")
        qq(gwas$P, main = cur_pheno)
        text(x=3,y=0.5,bquote(lambda==  .(lambda)))
        dev.off()
        # get rid of the data frame again
        rm(gwas,t,chisq,mchisq,lambda)
        count = count+1
      }
    }  
   }
  #merge the lambda list in one file per overall phenotype (thickness/surface / different covars) and save	
  final_lambda = do.call(rbind, lambdalist)
  write.csv(final_lambda,file=paste('bgenie_ukb43760_',pheno,'_',sample_set,'_lambdas_chisqStats.csv',sep=""),row.names = FALSE)
  lambdalist=list()
}	





