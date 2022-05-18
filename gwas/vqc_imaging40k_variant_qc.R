#-------------------------------------------------#
# Variant-level filtering for imaging data subset #
#-------------------------------------------------#

## This script is adapted by Dick Schijven from an original version written by Amaia Carrion Castillo
## Location of the original version: /data/workspaces/lag/workspaces/lg-ukbiobank/analysis/amaia/genetic_analysis/release_v3/imagingT1_N18057/snp_qc_subset_imaging_imp.R

options(stringsAsFactors = FALSE)

# Load required R packages
library(ggplot2)
library(grid)
library(gridExtra)
library(hexbin)

# Provide input arguments
args=commandArgs(trailingOnly=TRUE)
sourcefile=args[1] # Source file with configurations for the run
chr=args[2] # Chromosome number to process, provided through the command line.

source(sourcefile)

#-------------------------------------------#
# Define paths to relevant directories      #
#-------------------------------------------#

if (Sys.info()['sysname']=='Windows') {dir="P:/workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}

mfi_dir=paste(dir,"lg-ukbiobank/primary_data/genetic_data/snp/snp_release_v3/data/mfi", sep="")
hrc_dir=paste(dir,"lg-ukbiobank/derived_data/genetic_data/snp/snp_release_v2_qc/snp_QC/HRC",sep="")

# Make a QC directory within the subset directory (if it does not already exist)
if (dir.exists(working_dir)==FALSE) {
  dir.create(working_dir,showWarnings = FALSE)
}

setwd(working_dir)

# Define files to read
stats_f<-paste(stats_dir,list.files(stats_dir,pattern=paste(subset_name,"_chr",chr,".snpstats.txt",sep="")),sep="/") # SNPstats file
mfi_f<-paste(mfi_dir,list.files(mfi_dir,pattern=paste("chr",chr,"_v3*.*",sep="")),sep="/") # UKB SNPstats file (for total sample)
hrc_f<-paste(hrc_dir,list.files(hrc_dir,pattern=paste("chr",chr,".tab",sep="")),sep="/") # HRC statistics file


#-------------------------------------------#
# Read the data files and define columns    #
#-------------------------------------------#

# Imputation INFO and MAF from the UKB, release v3 (these statistics are calculated on the full genetic data)
mfi<-read.table(mfi_f, header = F)
# Set column names
colnames(mfi)<-c("SNP_ID","RS_ID","Position","A1.UKB","A2.UKB","MAF.UKB","MinorAllele.UKB","INFO.UKB")

# QCtool-based SNP statistics specific to the data subset
stats<-read.table(stats_f,header=T)
# Set column names
colnames(stats)<-c("SNP_ID", "RS_ID", "CHR", "Position", paste(colnames(stats)[5:ncol(stats)],".subset.QCtool",sep=""))
colnames(stats)[colnames(stats) == "NULL..subset.QCtool"] <- "NULL.subset.QCtool"

# HRC reference info (specific to the reference panel)
hrc<-read.table(hrc_f)
# Set column names
colnames(hrc)<-c("SNP_ID","RS_ID","AF.HRC")

# Split the information in the first column of the HRC file into multiple columns
hrc$CHR <- sapply(strsplit(hrc$SNP_ID,":"),"[[",1)
if( unique(hrc$CHR!="X")) {hrc$CHR <- as.numeric(hrc$CHR)}
hrc$Position <- as.numeric(sapply(strsplit(sapply(strsplit(hrc$SNP_ID,"_"),"[[",1),":"),"[[",2))
hrc$A1.HRC <- sapply(strsplit(hrc$SNP_ID,"_"),"[[",2)
hrc$A2.HRC <- sapply(strsplit(hrc$SNP_ID,"_"),"[[",3)
hrc$hrc<-1

# Create a MAF column for HRC, only 0-0.5
hrc$MAF.HRC <- hrc$AF.HRC
hrc$MAF.HRC[hrc$AF.HRC > 0.5] <- (1-hrc$AF.HRC[hrc$AF.HRC > 0.5])


#-------------------------------------------#
# Combine sources of information            #
#-------------------------------------------#

# Merge the UKB imputation information and the QCtool statistics
ukb <- merge(mfi,stats, by = c("RS_ID", "SNP_ID", "Position"), all=TRUE)

# Replace the SNP_ID column with a chr:pos:a1:a2 format column
ukb$SNP_ID <- paste(ukb$CHR,":",ukb$Position,"_",ukb$A1.UKB,"_",ukb$A2.UKB,sep="")


# for chr X, make hwe=female hwe
if (chr=="X") {
  ukb$HW_exact_p_value.subset.QCtool <- ukb$HW_females_exact_pvalue.subset.QCtool
}

# Combine with HRC SNP information
ukb_hrc <- merge(ukb, hrc, by=c("SNP_ID","CHR","Position"), all.x=TRUE, suffixes=c(".UKB",".HRC"))

# Check if HRC information for all available variants has been added to the dataframe
cat(paste("\n", sum(ukb_hrc$hrc == 1, na.rm = TRUE), " SNPs matched based on Chr:Position_Allele1_Allele2\n", sep=""))

# Sort by position
ukb_hrc <- ukb_hrc[with(ukb_hrc,order(CHR,Position)),]
ukb_hrc$hrc[is.na(ukb_hrc$hrc)] <- 0

# Print the number of SNPs (non-)matching HRC SNPs in the merged dataframe. Ideally the number of matching SNPs is equal to number of SNPs in the original HRC file.
cat(paste("\nNumber of SNPs matching HRC (1) in chromosome ",chr, sep=""))
table(ukb_hrc$hrc)
cat("\n")


#-------------------------------------------#
# Calculate differences between sources     #
#-------------------------------------------#

# MAF difference between the full UKB dataset and the HRC reference panel
ukb_hrc$MAFdiff_UKB_HRC <- ukb_hrc$MAF.UKB-ukb_hrc$MAF.HRC

# MAF difference between the HRC reference panel and the data subset
ukb_hrc$MAFdiff_HRC_subset <- ukb_hrc$MAF.HRC-ukb_hrc$minor_allele_frequency.subset.QCtool

# MAF difference between the full UKB dataset and the data subset
ukb_hrc$MAFdiff_UKB_subset <- ukb_hrc$MAF.UKB-ukb_hrc$minor_allele_frequency.subset.QCtool

# INFO difference between the full UKB dataset and the data subset
ukb_hrc$INFOdiff_UKB_subset <- ukb_hrc$INFO.UKB-ukb_hrc$info.subset.QCtool


#-------------------------------------------#
# Select and mark SNPs for keep/removal     #
#-------------------------------------------#

# First create a filter column and mark all SNPs with 0 (means all SNPs will be kept)
ukb_hrc$filter <- 0

## Evaluate all filtering criteria and mark variants with 1 (removal) if selected:

cat("Applying filtering:\n")

# MAF threshold
if (!(is.null(maf_thr))){
  w <- which(ukb_hrc$minor_allele_frequency.subset.QCtool < maf_thr)
  ukb_hrc$filter[w] <- 1
  cat(paste("Variants with minor_allele_frequency.subset.QCtool < ", maf_thr, ": ", length(w), "\n", sep=""))
  rm(w)
}

# MAF difference
if (!(is.null(maf_diff))){
  w <- which(ukb_hrc$MAFdiff_UKB_subset > maf_diff | ukb_hrc$MAFdiff_UKB_subset < -maf_diff | ukb_hrc$MAFdiff_HRC_subset > maf_diff | ukb_hrc$MAFdiff_HRC_subset < -maf_diff)
  ukb_hrc$filter[w] <- 1
  cat(paste("Variants with MAFdiff_UKB_subset < ", -maf_diff, " or > ", maf_diff, ", and with MAFdiff_HRC_subset < ", -maf_diff, " or > ", maf_diff, ": ", length(w), "\n", sep=""))
  rm(w)
}

# INFO threshold
if (!(is.null(info_thr))){
  w <- which(ukb_hrc$info.subset.QCtool < info_thr)
  ukb_hrc$filter[w] <- 1
  cat(paste("Variants with info.subset.QCtool < ", info_thr, ": ", length(w), "\n", sep=""))
  rm(w)
}

#INFO difference
if (!(is.null(info_diff))){
  w <- which(ukb_hrc$INFOdiff_UKB_subset < -info_diff | ukb_hrc$INFOdiff_UKB_subset > info_diff)
  ukb_hrc$filter[w] <- 1
  cat(paste("Variants with INFOdiff_UKB_subset < ", -info_diff, " or > ", info_diff, ": ", length(w), "\n", sep=""))
  rm(w)
}

# HWE threshold
if (!(is.null(hwe_thr))){
  w <- which(ukb_hrc$HW_exact_p_value.subset.QCtool < hwe_thr)
  ukb_hrc$filter[w] <- 1
  cat(paste("Variants with HW_exact_p_value.subset.QCtool < ", hwe_thr, ": ", length(w), "\n", sep=""))
  rm(w)
}

# Variant-level missingness
if (!(is.null(geno_miss))){
  w <- which(ukb_hrc$missing_proportion.subset.QCtool > geno_miss)
  ukb_hrc$filter[w] <- 1
  cat(paste("Variants with missing_proportion.subset.QCtool > ", geno_miss, ": ", length(w), "\n", sep=""))
  rm(w)
}

# Remove indels
if (remove_indel == T & !(is.null(remove_indel))){
  w <- which( (ukb_hrc$A1.UKB != "A" & ukb_hrc$A1.UKB != "T" & ukb_hrc$A1.UKB != "C" & ukb_hrc$A1.UKB != "G") | (ukb_hrc$A2.UKB != "A" & ukb_hrc$A2.UKB != "T" & ukb_hrc$A2.UKB != "C" & ukb_hrc$A2.UKB != "G") )
  ukb_hrc$filter[w] <- 1
  cat(paste("Variants other than SNPs (such as indels): ", length(w), "\n", sep=""))
  rm(w)
}

# Remove multiallelic variants
if (remove_multiallelic == T & !(is.null(remove_multiallelic))){
  n_occur <- data.frame(table(ukb_hrc$RS_ID.UKB))
  duplicates <- n_occur$Var1[n_occur$Freq > 1]
  w <- which(ukb_hrc$RS_ID.UKB %in% duplicates)
  ukb_hrc$filter[w] <- 1
  cat(paste("Multiallelic variants: ", length(duplicates), " unique variants, ", length(w), " major-minor allele pairs", "\n", sep=""))
  rm(w)
}

# Show an overview of the numbers of SNPs to keep and remove
cat("\nSummary of variants to keep (0) and to exclude (1)")
table(ukb_hrc$filter)
cat("\n")


#---------------------------------------------#
# Write tables of variants to keep and remove #
#---------------------------------------------#

if (compact_output == FALSE){
  
  # Save a full output file with all SNP statistics from different sources
  out_f <- gsub(stats_dir,working_dir,gsub(".txt","_mfi_hrc.txt",stats_f))
  
  cols2remove <- c("comment.subset.QCtool","impute_info.subset.QCtool")
  cols2keep <- colnames(ukb_hrc)[!(colnames(ukb_hrc) %in% cols2remove)]
  
  write.table(subset(ukb_hrc, select=cols2keep),out_f,sep="\t",row.names=FALSE,quote=FALSE)
  
  rm(cols2keep, out_f)
  
  # Save detailed lists of variants to keep and remove
  cols2keep<-c("SNP_ID","CHR","Position","RS_ID.UKB","A1.UKB","A2.UKB", "minor_allele.subset.QCtool", "minor_allele_frequency.subset.QCtool",
               "info.subset.QCtool", "HW_exact_p_value.subset.QCtool", "missing_proportion.subset.QCtool", "hrc", "MAFdiff_HRC_subset",
               "MAFdiff_UKB_subset", "INFOdiff_UKB_subset")
  
  keep_f <- gsub(stats_dir, working_dir, gsub(".txt","_mfi_hrc.snps2keep",stats_f))
  remove_f <- gsub(stats_dir, working_dir, gsub(".txt","_mfi_hrc.snps2remove",stats_f))
  write.table(subset(ukb_hrc, filter==0, select=cols2keep), keep_f, sep="\t", row.names=FALSE, quote=FALSE)
  write.table(subset(ukb_hrc, filter==1, select=cols2keep), remove_f, sep="\t", row.names=FALSE, quote=FALSE)
  
  rm(cols2keep, keep_f, remove_f)
  
}

if (compact_output == TRUE) {
  
  # Only save a list of variant IDs to keep and remove
  cols2keep<-c("SNP_ID", "RS_ID.UKB")
  
  keep_f <- gsub(stats_dir, working_dir, gsub(".txt","_mfi_hrc.compact.snps2keep",stats_f))
  remove_f <- gsub(stats_dir, working_dir, gsub(".txt","_mfi_hrc.compact.snps2remove",stats_f))
  write.table(subset(ukb_hrc, filter==0, select=cols2keep), keep_f, sep="\t", col.names = TRUE, row.names=FALSE, quote=FALSE)
  write.table(subset(ukb_hrc, filter==1, select=cols2keep), remove_f, sep="\t", col.names = TRUE, row.names=FALSE, quote=FALSE)
  
  rm(cols2keep, keep_f, remove_f)
  
}


#-------------------------------------------#
# Plots (if plots == T)                     #
#-------------------------------------------#

if (plots == T){
  
  # Within the QC directory, create a plots directory (if it does not already exist)/
  if (dir.exists(paste(working_dir,"plots",sep="/"))==FALSE) {
    dir.create(path=paste(working_dir,"plots",sep="/"),showWarnings = FALSE)
  }

  # Make the first plot: hexbins of info scores in the full UKB against the subset
  INFO_plot1<-ggplot(ukb_hrc,aes(x=INFO.UKB,y=impute_info.subset.QCtool)) + 
    geom_hex(bins=50) +
    scale_fill_gradient(trans = "log10") +
    theme_bw() + 
    geom_vline(xintercept=info_thr,color="red") +
    geom_hline(yintercept=info_thr,color="red") +
    ggtitle(paste("INFO comparison: chr ",chr,"\nTotal UKB vs UKB subset\nUKB subset: ",subset_name,sep="")) +
    xlab("INFO - UKB Total") +
    ylab("INFO - UKB Subset")
  
  # Save the first plot
  plot_f<-paste(working_dir,"/plots/INFO_UKBtotal_UKBsubset_",subset_name,"_chr",chr,".png",sep="")
  png(plot_f)
  print(INFO_plot1)
  dev.off()
  
  rm(plot_f, INFO_plot1)
  
  # Make the second plot: distribution of INFO differences between the full UKB and the subset
  INFO_plot2<-ggplot(ukb_hrc,aes(x=INFOdiff_UKB_subset)) + 
    geom_histogram() +
    theme_bw() +
    ggtitle(paste("INFO difference: chr ",chr,"\nTotal UKB vs UKB subset\nUKB subset: ",subset_name,sep="")) +
    xlab("INFO difference") +
    ylab("Count")
  
  # Save the second plot
  plot_f2<-paste(working_dir,"/plots/INFOdiff_UKBtotal_UKBsubset_",subset_name,"_chr",chr,".png",sep="")
  png(plot_f2)
  print(INFO_plot2)
  dev.off()
  
  rm(plot_f2,INFO_plot2)
  
  ## Compare Allele frequencies: total UKB, subset UKB and HRC
  # Make the first plot
  AF_plot1 <- ggplot(ukb_hrc,aes(x=MAF.UKB,y=minor_allele_frequency.subset.QCtool)) +
    geom_point() +
    theme_bw() +
    ggtitle(paste("UKB Subset vs. UKB Total\nUKB Subset: ",subset_name,sep="")) +
    xlab("MAF - UKB Total") +
    ylab("MAF - UKB Subset")
  
  # Make the second plot
  AF_plot2 <- ggplot(ukb_hrc,aes(x=MAF.UKB,y=MAF.HRC)) +
    geom_point() +
    theme_bw() +
    ggtitle(paste("HRC vs. Total UKB",sep="")) + 
    xlab("MAF - UKB Total") +
    ylab("MAF - HRC")
  
  # Make the third plot
  AF_plot3 <- ggplot(ukb_hrc,aes(x=minor_allele_frequency.subset.QCtool,y=MAF.HRC)) +
    geom_point() +
    theme_bw() +
    ggtitle(paste("HRC vs UKB Subset\nUKB Subset: ",subset_name,sep="")) +
    xlab("MAF - UKB Subset") +
    ylab("MAF - HRC")
  
  # Combine the three plots and save
  plot_f<-paste(working_dir,"/plots/AF_HRC_UKBtotal_UKBsubset_",subset_name,"_chr",chr,".png",sep="")
  png(plot_f)
  grid.arrange(AF_plot1,AF_plot2,AF_plot3, top=paste("Allele frequency comparison: chr ", chr,"\nNumber of SNPs: ",NROW(ukb_hrc),sep=""))
  dev.off()
  
  rm(AF_plot1,AF_plot2,AF_plot3,plot_f)
  
  
  
  ## PLOTS WITH ONLY VARIANTS THAT ARE MARKED TO KEEP ##
  
  ## Compare INFO fields: UKB MFI file vs QCtool INFO for the subset
  
  # Make the first plot: hexbins of info scores in the full UKB against the subset
  INFO_plot1<-ggplot(subset(ukb_hrc,filter==0),aes(x=INFO.UKB,y=impute_info.subset.QCtool)) + 
    geom_hex(bins=50) +
    scale_fill_gradient(trans = "log10") +
    theme_bw() + 
    geom_vline(xintercept=info_thr,color="red") +
    geom_hline(yintercept=info_thr,color="red") +
    ggtitle(paste("INFO comparison: chr ",chr,"\nTotal UKB vs UKB subset\nUKB subset: ",subset_name,sep="")) +
    xlab("INFO - UKB Total") +
    ylab("INFO - UKB Subset")
  
  # Save the first plot
  plot_f<-paste(working_dir,"/plots/INFO_UKBtotal_UKBsubset_",subset_name,"_chr",chr,"_snps2keep.png",sep="")
  png(plot_f)
  print(INFO_plot1)
  dev.off()
  
  rm(plot_f,INFO_plot1)
  
  # Make the second plot: distribution of INFO differences between the full UKB and the subset
  INFO_plot2<-ggplot(subset(ukb_hrc,filter==0),aes(x=INFOdiff_UKB_subset)) + 
    geom_histogram() +
    theme_bw() +
    ggtitle(paste("INFO difference: chr ",chr,"\nTotal UKB vs UKB subset\nUKB subset: ",subset_name,sep="")) +
    xlab("INFO difference") +
    ylab("Count")
  
  # Save the second plot
  plot_f2<-paste(working_dir,"/plots/INFOdiff_UKBtotal_UKBsubset_",subset_name,"_chr",chr,"_snps2keep.png",sep="")
  png(plot_f2)
  print(INFO_plot2)
  dev.off()
  
  rm(plot_f2,INFO_plot2)
  
  ## Compare Allele frequencies: total UKB, subset UKB and HRC
  # Make the first plot
  AF_plot1 <- ggplot(subset(ukb_hrc,filter==0),aes(x=MAF.UKB,y=minor_allele_frequency.subset.QCtool)) +
    geom_point() +
    theme_bw() +
    ggtitle(paste("UKB Subset vs. UKB Total\nUKB Subset: ",subset_name,sep="")) +
    xlab("MAF - UKB Total") +
    ylab("MAF - UKB Subset")
  
  # Make the second plot
  AF_plot2 <- ggplot(subset(ukb_hrc,filter==0),aes(x=MAF.UKB,y=MAF.HRC)) +
    geom_point() +
    theme_bw() +
    ggtitle(paste("HRC vs. Total UKB",sep="")) + 
    xlab("MAF - UKB Total") +
    ylab("MAF - HRC")
  
  # Make the third plot
  AF_plot3 <- ggplot(subset(ukb_hrc,filter==0),aes(x=minor_allele_frequency.subset.QCtool,y=MAF.HRC)) +
    geom_point() +
    theme_bw() +
    ggtitle(paste("HRC vs UKB Subset\nUKB Subset: ",subset_name,sep="")) +
    xlab("MAF - UKB Subset") +
    ylab("MAF - HRC")
  
  # Combine the three plots and save
  plot_f <- paste(working_dir,"/plots/AF_HRC_UKBtotal_UKBsubset_",subset_name,"_chr",chr,"_snps2keep.png",sep="")
  png(plot_f)
  grid.arrange(AF_plot1,AF_plot2,AF_plot3, top = paste("Allele frequency comparison: chr ", chr,"\nNumber of SNPs: ", nrow(subset(ukb_hrc,filter==0)), sep=""))
  dev.off()
  
  rm(AF_plot1,AF_plot2,AF_plot3,plot_f)
}
