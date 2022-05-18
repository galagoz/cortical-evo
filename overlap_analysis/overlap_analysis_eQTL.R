#' ---
#' title: "overlap_analysis_eQTL"
#' author: "Gokberk Alagoz"
#' date: "September 26, 2021"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=FALSE)
library(GenomicRanges)
library(biomaRt)

#' 
#' ## R Markdown
#' 
## ----paths, echo=FALSE---------------------------------------------------------------------------------------------------------------------------
bedfile = args[1] # "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/GRCh37/neanDepRegions_hg19.sorted.bed"
clumpedSumstatsDir = args[2] #"/data/clusterfs/lag/users/gokala/enigma-evol/eqtl/clumped_sumstats/european_lr"
outDir = args[3] # "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/eqtl/european_lr/results"
genotypeF = "/data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink/1KG_phase3_GRCh37_EUR_nonFIN_allchr"

#' 
#' ## Including Plots
#' 
## ----parse, echo=FALSE---------------------------------------------------------------------------------------------------------------------------

# eQTL data downloaded from PsychENCODE
feqtl = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/old_results/eqtl/DER-08a_hg19_eQTL.significant.txt"

# Parse to get trait name
clumpfileloc = dir(clumpedSumstatsDir, pattern = ".clumped", full.names = T)

phenoname = paste0(sapply(clumpfileloc, function (x) {unlist(strsplit(x, "_", fixed = TRUE))[5]}), "_", 
                   sapply(clumpfileloc, function (x) {unlist(strsplit(x, "_", fixed = TRUE))[6]}), "_", 
                   sapply(clumpfileloc, function (x) {unlist(strsplit(x, "_", fixed = TRUE))[7]}))

annot_name = unlist(strsplit(unlist(strsplit(bedfile, "/", fixed=T))[15], ".", fixed=T))[1]

#' 
## ----fullSA_left---------------------------------------------------------------------------------------------------------------------------------

# Full surface area (LEFT HEM.)
# Start by getting all the clumped SNPs only for full surface area

fullsurfind = which(phenoname=="surface_Full_le")
clump = read.table(clumpfileloc[fullsurfind],header=TRUE)

# Loop over all clumped SNPs
#i=1
for (i in 1:nrow(clump)) {
    
    # Find all SNPs in LD (r2>0.6) with each clumped SNP
  
    system(paste0("module load plink/1.9b6 \
                   plink --bfile ", genotypeF, " --r2 --ld-window-kb 10000 --ld-window 2000 --ld-window-r2 0.6 --ld-snp ",
                  clump$SNP[i], " --out tmpld_", annot_name))
    
    # Read in the LD calculated from plink
    
    LD = read.table(paste0("tmpld_",annot_name,".ld"),header=TRUE)
    
    # Turn into a genomic ranges object
    
    if (i==1) { 
       LDSNPs = GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP = LD$SNP_B,indexSNP = LD$SNP_A)
    } else {
       LDSNPs = c(GRanges(LD$CHR_B,IRanges(LD$BP_B,LD$BP_B),SNP = LD$SNP_B,indexSNP = LD$SNP_A), LDSNPs)
    }
}
  
# Read in the eQTL data
eqtl = read.table(feqtl, header=TRUE)
eqtl.GR = GRanges(gsub("chr", "", eqtl$SNP_chr), IRanges(eqtl$SNP_start, eqtl$SNP_end))
mcols(eqtl.GR) = eqtl[, c(1:8, 12:15)]

# Read in the BED file containing the annotation file
annot = read.table(bedfile, header=FALSE)
annotGR = GRanges(gsub("chr", "", annot$V1), IRanges(annot$V2, annot$V3), type = annot$V4)
        
# Find overlaps
olap = findOverlaps(annotGR, LDSNPs)
regionalSA_vars = unique(LDSNPs$indexSNP[subjectHits(olap)])
regionalAnnot = LDSNPs[which(!is.na(match(LDSNPs$indexSNP, regionalSA_vars)))]

# Outputting the number of loci that overlap with an HGE
cat('Number of loci that overlap with ', annot_name, ' is: ', length(unique(regionalAnnot$indexSNP)), '\n')

olap2 = findOverlaps(eqtl.GR, regionalAnnot)
regionalSA_eqtl_vars  = unique(regionalAnnot$indexSNP[subjectHits(olap2)])
cat('Number of ', annot_name, ' overlapping loci that also have an eQTL: ', length(unique(regionalSA_eqtl_vars$indexSNP[subjectHits(olap2)])), '\n')
eqtlgenes = unique(eqtl.GR$gene_id[queryHits(olap2)])
        
# Remove the . annotation
eqtlgenes = sapply(eqtlgenes, function (x) {unlist(strsplit(x, ".", fixed=TRUE))[1]})

if (length(unique(regionalSA_eqtl_vars$indexSNP[subjectHits(olap2)])) > 0) {
  
  # Convert these genes to hgnc_id
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                 dataset="hsapiens_gene_ensembl",
                 host="feb2014.archive.ensembl.org")
  geneannot = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                    filters = "ensembl_gene_id",
                    values = eqtlgenes,
                    mart = mart)
  cat('These eQTLs impact ', length(which(geneannot == "protein_coding")), ' protein-coding eGenes\n')
  write.csv(geneannot, file = paste0(outDir, "/leftHem_SA_", annot_name, "_lr_eqtl_genes.csv"),
            row.names = FALSE, quote = FALSE)
  
} else {
  print("There are not any CSA-associated variant - eQTL overlap")
}

if (length(regionalSA_vars) > 0) {

  m = matrix(NA, nrow = length(regionalSA_vars), ncol = 3)
  d = as.data.frame(m)
  colnames(d) = c("olap_snps", "eqtl_snps", "region")
  
  d$olap_snps = regionalSA_vars
  d$eqtl_snps[match(regionalSA_eqtl_vars, d$olap_snps)] = "Yes"
  d$region = "full_leftHem"
  d$olap_snps[match(regionalSA_eqtl_vars, d$olap_snps)]
  write.csv(d, file = paste0(outDir, "/leftHem_SA_", annot_name, "_lr_olap_snps.csv"),
            row.names = FALSE, quote = FALSE)
} else {
  print("There are not any CSA-associated variants - annotation overlap")
}


#' 
## ----fullSA_right--------------------------------------------------------------------------------------------------------------------------------

# Full surface area (RIGHT HEM.)
# Start by getting all the clumped SNPs only for full surface area

fullsurfind = which(phenoname == "surface_re_Full")
clump = read.table(clumpfileloc[fullsurfind], header = TRUE)

# Loop over all clumped SNPs

for (i in 1:nrow(clump)) {
  
  # Find all SNPs in LD (r2>0.6) with each clumped SNP
  
    system(paste0("module load plink/1.9b6 \
                   plink --bfile ", genotypeF, " --r2 --ld-window-kb 10000 --ld-window 2000 --ld-window-r2 0.6 --ld-snp ",
                  clump$SNP[i], " --out tmpld_", annot_name))
  
  # Read in the LD calculated from plink
  
  LD = read.table(paste0("tmpld_", annot_name, ".ld"), header = TRUE)
  
  # Turn into a genomic ranges object
  
  if (i==1) { 
    LDSNPs = GRanges(LD$CHR_B, IRanges(LD$BP_B, LD$BP_B), SNP = LD$SNP_B, indexSNP = LD$SNP_A)
  } else {
    LDSNPs = c(GRanges(LD$CHR_B, IRanges(LD$BP_B, LD$BP_B), SNP = LD$SNP_B, indexSNP = LD$SNP_A), LDSNPs)
  }
}

# Read in the eQTL data
eqtl = read.table(feqtl, header = TRUE)
eqtl.GR = GRanges(gsub("chr", "", eqtl$SNP_chr), IRanges(eqtl$SNP_start, eqtl$SNP_end))
mcols(eqtl.GR) = eqtl[, c(1:8,12:15)]

# Read in the BED file containing the annotation file
annot = read.table(bedfile, header = FALSE)
annotGR = GRanges(gsub("chr", "", annot$V1), IRanges(annot$V2, annot$V3), type = annot$V4)

# Find overlaps
olap = findOverlaps(annotGR, LDSNPs)
regionalSA_vars = unique(LDSNPs$indexSNP[subjectHits(olap)])
regionalAnnot = LDSNPs[which(!is.na(match(LDSNPs$indexSNP, regionalSA_vars)))]

# Outputting the number of loci that overlap with an HGE
cat('Number of loci that overlap with ', annot_name, ' is: ', length(unique(regionalAnnot$indexSNP)), '\n')

olap2 = findOverlaps(eqtl.GR, regionalAnnot)
regionalSA_eqtl_vars = unique(regionalAnnot$indexSNP[subjectHits(olap2)])
cat('Number of ', annot_name, ' overlapping loci that also have an eQTL: ', length(unique(regionalSA_eqtl_vars$indexSNP[subjectHits(olap2)])), '\n')
eqtlgenes = unique(eqtl.GR$gene_id[queryHits(olap2)])

# Remove the . annotation
eqtlgenes = sapply(eqtlgenes, function (x) {unlist(strsplit(x, ".", fixed = TRUE))[1]})

if (length(unique(regionalSA_eqtl_vars$indexSNP[subjectHits(olap2)])) > 0) {
  
  # Convert these genes to hgnc_id
  mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                 dataset = "hsapiens_gene_ensembl",
                 host = "feb2014.archive.ensembl.org")
  geneannot = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                    filters = "ensembl_gene_id",
                    values = eqtlgenes,
                    mart = mart)
  cat('These eQTLs impact ', length(which(geneannot == "protein_coding")), ' protein-coding eGenes\n')
  write.csv(geneannot, file = paste0(outDir, "/rightHem_SA_", annot_name, "_lr_eqtl_genes.csv"), row.names = FALSE, quote = FALSE)
  
} else {
  print("There are not any CSA-associated variant - eQTL overlap")
}

if (length(regionalSA_vars) > 0) {

  m = matrix(NA, nrow = length(regionalSA_vars), ncol = 3)
  d = as.data.frame(m)
  colnames(d) = c("olap_snps", "eqtl_snps", "region")
  
  d$olap_snps = regionalSA_vars
  d$eqtl_snps[match(regionalSA_eqtl_vars, d$olap_snps)] = "Yes"
  d$region = "full_rightHem"
  d$olap_snps[match(regionalSA_eqtl_vars, d$olap_snps)]
  write.csv(d, file = paste0(outDir, "/rightHem_SA_", annot_name, "_lr_olap_snps.csv"),
            row.names = FALSE, quote = FALSE)
} else {
  print("There are not any CSA-associated variants - annotation overlap")
}


#' 
## ----regionalSA_left-----------------------------------------------------------------------------------------------------------------------------

# All regional surface areas (LEFT HEM.)
# Start by getting all the clumped SNPs only for any regional surface area

surfareaind = grep("surface_le", phenoname) 
surfareaind = surfareaind[surfareaind != which(phenoname == "surface_le_Full")]

# Initiate a clump file, we will add other clump files to this one
clump = read.table(clumpfileloc[surfareaind[1]], header = TRUE)

first_region = strsplit(clumpfileloc[surfareaind[1]],"_",fixed=TRUE)[[1]][8]
clump$region = first_region

for (j in 2:length(surfareaind)) {
    
  if (file.exists(clumpfileloc[surfareaind[j]])) {
    tmp_clump = read.table(clumpfileloc[surfareaind[j]], header=TRUE)
    tmp_region_name = strsplit(clumpfileloc[surfareaind[j]],"_",fixed=TRUE)[[1]][8]
    tmp_clump$region = tmp_region_name
    clump = rbind(clump, tmp_clump, header=TRUE)
  }
}

clump = clump[clump$SNP!=TRUE,]

  # Loop over all SNPs
  for (i in 1:nrow(clump)) {
    
    # Find all SNPs in LD (r2>0.6) with each clumped SNP
    system(paste0("module load plink/1.9b6 \
                   plink --bfile ", genotypeF, " --r2 --ld-window-kb 10000 --ld-window 2000 --ld-window-r2 0.6 --ld-snp ",
                  clump$SNP[i], " --out tmpld_", annot_name))
    
    # Read in the LD calculated from plink
    LD = read.table(paste0("reg_tmpld_", annot_name, ".ld"), header = TRUE)
    
    # Turn into a genomic ranges object
    if (i==1) { 
      LDSNPs = GRanges(LD$CHR_B, IRanges(LD$BP_B, LD$BP_B), SNP=LD$SNP_B, indexSNP=LD$SNP_A)
    } else {
      LDSNPs = c(GRanges(LD$CHR_B, IRanges(LD$BP_B, LD$BP_B), SNP=LD$SNP_B, indexSNP=LD$SNP_A), LDSNPs)
    }
  }
  
# Make sure that all index SNPs are in this list
LDSNPs = c(GRanges(clump$CHR, IRanges(clump$BP, clump$BP), SNP = clump$SNP, indexSNP = clump$SNP), LDSNPs)
  
#save(LDSNPs,file="RegionalLDSNPs.Rdata")

# Find overlaps
olap = findOverlaps(annotGR, LDSNPs)
regionalSA_vars = unique(LDSNPs$indexSNP[subjectHits(olap)])
regionalAnnot = LDSNPs[which(!is.na(match(LDSNPs$indexSNP, regionalSA_vars)))]

# Outputting the number of loci that overlap with an HGE
cat('Number of loci that overlap with ', annot_name, ' is: ', length(unique(regionalAnnot$indexSNP)), '\n')

olap2 = findOverlaps(eqtl.GR, regionalAnnot)
regionalSA_eqtl_vars  =unique(regionalAnnot$indexSNP[subjectHits(olap2)])
cat('Number of ', annot_name, ' overlapping loci that also have an eQTL: ', length(unique(regionalAnnot$indexSNP[subjectHits(olap2)])), '\n')
eqtlgenes = unique(eqtl.GR$gene_id[queryHits(olap2)])
  
# Remove the . annotation
eqtlgenes = sapply(eqtlgenes, function (x) {unlist(strsplit(x, ".", fixed = TRUE))[1]})

if (length(unique(regionalAnnot$indexSNP[subjectHits(olap2)])) > 0) {
  
  # Convert these genes to hgnc_id
  mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                 dataset = "hsapiens_gene_ensembl",
                 host = "feb2014.archive.ensembl.org")
  geneannot = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                    filters = "ensembl_gene_id", values = eqtlgenes,mart = mart)
  cat('These eQTLs impact ', length(which(geneannot == "protein_coding")), ' protein-coding eGenes\n')
  write.csv(geneannot, file = paste0(outDir, "/regionalSA_le_", annot_name, "_lr_eqtl_genes.csv"), row.names = FALSE, quote = FALSE)

} else {
  print("There are not any CSA-associated variant - eQTL overlap")
}

if (length(regionalSA_vars) > 0) {

  m = matrix(NA, nrow = length(regionalSA_vars), ncol = 3)
  d = as.data.frame(m)
  colnames(d) = c("olap_snps", "eqtl_snps", "region")
  
  d$olap_snps = regionalSA_vars
  d$eqtl_snps[match(regionalSA_eqtl_vars, d$olap_snps)] = "Yes"
  d$region = clump[match(regionalSA_vars, clump$SNP),]$region
  
  write.csv(d, file = paste0(outDir, "/regionalSA_le_", annot_name, "_lr_olap_snps.csv"),
            row.names = FALSE, quote = FALSE)
} else {
  print("There are not any CSA-associated variants - annotation overlap")
}


#' 
## ----regionalSA_right----------------------------------------------------------------------------------------------------------------------------

# All regional surface areas (RIGHT HEM.)
# Start by getting all the clumped SNPs only for any regional surface area

surfareaind = grep("surface_re", phenoname)
surfareaind = surfareaind[surfareaind != which(phenoname == "surface_re_Full")]

# Initiate a clump file, we will add other clump files to this one
clump = read.table(clumpfileloc[surfareaind[1]], header = TRUE)

first_region = strsplit(clumpfileloc[surfareaind[1]],"_",fixed=TRUE)[[1]][8]
clump$region = first_region

for (j in 2:length(surfareaind)) {
    
  if (file.exists(clumpfileloc[surfareaind[j]])) {
    tmp_clump = read.table(clumpfileloc[surfareaind[j]], header=TRUE)
    tmp_region_name = strsplit(clumpfileloc[surfareaind[j]],"_",fixed=TRUE)[[1]][8]
    tmp_clump$region = tmp_region_name
    clump = rbind(clump, tmp_clump, header=TRUE)
  }
}

clump = clump[clump$SNP!=TRUE,]

# Loop over all SNPs
for (i in 1:nrow(clump)) {
  
  # Find all SNPs in LD (r2>0.6) with each clumped SNP
    system(paste0("module load plink/1.9b6 \
                   plink --bfile ", genotypeF, " --r2 --ld-window-kb 10000 --ld-window 2000 --ld-window-r2 0.6 --ld-snp ",
                  clump$SNP[i], " --out tmpld_", annot_name))
  
  # Read in the LD calculated from plink
  LD = read.table(paste0("reg_tmpld_", annot_name, ".ld"), header = TRUE)
  
  # Turn into a genomic ranges object
  if (i==1) { 
    LDSNPs = GRanges(LD$CHR_B, IRanges(LD$BP_B, LD$BP_B), SNP = LD$SNP_B, indexSNP = LD$SNP_A)
  } else {
    LDSNPs = c(GRanges(LD$CHR_B, IRanges(LD$BP_B, LD$BP_B), SNP = LD$SNP_B, indexSNP = LD$SNP_A), LDSNPs)
  }
}

# Make sure that all index SNPs are in this list
LDSNPs = c(GRanges(clump$CHR, IRanges(clump$BP, clump$BP), SNP = clump$SNP, indexSNP = clump$SNP), LDSNPs)

#save(LDSNPs,file="RegionalLDSNPs.Rdata")

# Find overlaps
olap = findOverlaps(annotGR, LDSNPs)
regionalSA_vars = unique(LDSNPs$indexSNP[subjectHits(olap)])
regionalAnnot = LDSNPs[which(!is.na(match(LDSNPs$indexSNP, regionalSA_vars)))]

# Outputting the number of loci that overlap with an HGE
cat('Number of loci that overlap with ', annot_name, ' is: ', length(unique(regionalAnnot$indexSNP)), '\n')

olap2 = findOverlaps(eqtl.GR, regionalAnnot)
regionalSA_eqtl_vars  =unique(regionalAnnot$indexSNP[subjectHits(olap2)])
cat('Number of ', annot_name, ' overlapping loci that also have an eQTL: ', length(unique(regionalAnnot$indexSNP[subjectHits(olap2)])), '\n')
eqtlgenes = unique(eqtl.GR$gene_id[queryHits(olap2)])

# Remove the . annotation
eqtlgenes = sapply(eqtlgenes, function (x) {unlist(strsplit(x, ".", fixed = TRUE))[1]})

if (length(unique(regionalAnnot$indexSNP[subjectHits(olap2)])) > 0) {
  
  # Convert these genes to hgnc_id
  mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                 dataset = "hsapiens_gene_ensembl",
                 host = "feb2014.archive.ensembl.org")
  geneannot = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                    filters = "ensembl_gene_id", values = eqtlgenes, mart = mart)
  cat('These eQTLs impact ', length(which(geneannot == "protein_coding")), ' protein-coding eGenes\n')
  write.csv(geneannot, file = paste0(outDir, "/regionalSA_ri_", annot_name, "_lr_eqtl_genes.csv"), row.names = FALSE, quote = FALSE)
  
} else {
  print("There are not any CSA-associated variant - eQTL overlap")
}

if (length(regionalSA_vars) > 0) {

  m = matrix(NA, nrow = length(regionalSA_vars), ncol = 3)
  d = as.data.frame(m)
  colnames(d) = c("olap_snps", "eqtl_snps", "region")
  
  d$olap_snps = regionalSA_vars
  d$eqtl_snps[match(regionalSA_eqtl_vars, d$olap_snps)] = "Yes"
  d$region = clump[match(regionalSA_vars, clump$SNP),]$region
  
  write.csv(d, file = paste0(outDir, "/regionalSA_ri_", annot_name, "_lr_olap_snps.csv"),
            row.names = FALSE, quote = FALSE)
} else {
  print("There are not any CSA-associated variants - annotation overlap")
}


#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()
#knitr::purl("overlap_analysis_eQTL.Rmd", "overlap_analysis_eQTL.R", documentation = 2) # save .Rmd as an .R file to submit to Grid

