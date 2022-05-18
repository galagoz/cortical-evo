# ========================================================
# Make tables of LDSC partitioned heritability results

# The results files (.results) from the LDSC run should be organized
# by annotation, with separate directories for each annotation. 

# Updated for the nonGC version

# ========================================================

library(tidyverse)
library(openxlsx)
options(stringsAsFactors=FALSE)

annots = dir(path = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/final_results/european_lr/results", full.names = F, recursive = F, pattern = "sorted")

#####
# Define a custom p.adjust to correct for the total number of independent DTI tracts.

my.p.adj = function (p, method = p.adjust.methods, n = length(p)) 
{
  method <- match.arg(method)
  if (method == "fdr") 
    method <- "BH"
  nm <- names(p)
  p <- as.numeric(p)
  p0 <- setNames(p, nm)
  if (all(nna <- !is.na(p))) 
    nna <- TRUE
  p <- p[nna]
  lp <- length(p)
  #stopifnot(n >= lp)
  if (n <= 1) 
    return(p0)
  if (n == 2 && method == "hommel") 
    method <- "hochberg"
  p0[nna] <- switch(method, bonferroni = pmin(1, n * p), holm = {
    i <- seq_len(lp)
    o <- order(p)
    ro <- order(o)
    pmin(1, cummax((n + 1L - i) * p[o]))[ro]
  }, hommel = {
    if (n > lp) p <- c(p, rep.int(1, n - lp))
    i <- seq_len(n)
    o <- order(p)
    p <- p[o]
    ro <- order(o)
    q <- pa <- rep.int(min(n * p/i), n)
    for (j in (n - 1L):2L) {
      ij <- seq_len(n - j + 1L)
      i2 <- (n - j + 2L):n
      q1 <- min(j * p[i2]/(2L:j))
      q[ij] <- pmin(j * p[ij], q1)
      q[i2] <- q[n - j + 1L]
      pa <- pmax(pa, q)
    }
    pmax(pa, p)[if (lp < n) ro[1L:lp] else ro]
  }, hochberg = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin((n + 1L - i) * p[o]))[ro]
  }, BH = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin(n/i * p[o]))[ro]
  }, BY = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    q <- sum(1/(1L:n))
    pmin(1, cummin(q * n/i * p[o]))[ro]
  }, none = p)
  p0
}

#####

#i=2
for (i in 1:length(annots)){
      print(annots[i])
      files = Sys.glob(path = paste0("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/final_results/european_lr/results/",annots[i],"/*.gz.results"))
      partheritresults = data.frame(Category = character(0),
                                    Prop._SNPs= numeric(0),
                                    Prop._h2= numeric(0),
                                    Prop._h2_std_error= numeric(0),
                                    Enrichment= numeric(0),
                                    Enrichment_std_error= numeric(0),
                                    Enrichment_p= numeric(0),
                                    Annotation=character(0),
                                    Analysis=character(0),
                                    Region=character(0)) #Will have a matrix with rows = number of E3MAs and columns = # of annotations
      #j=7
      for (j in 1:length(files)) {
        results = read.table(files[j],header=TRUE)
        info1 = str_split(files[j], pattern = "/")
        info2 = str_split(info1[[1]][17], pattern = "_")
        results$Annotation = info1[[1]][16]
        results$Analysis = if_else(grepl("hick",info2[[1]][4]), "Thickness", "Surface Area")
        #results$Analysis = "WMtracts"
        results$Region = paste(info2[[1]][5:6],collapse = "_")
        partheritresults = rbind(partheritresults,results[1,])
      }
      partheritresults = partheritresults %>% 
        group_by(Analysis) %>%
        mutate(fdr = my.p.adj(Enrichment_p, method = "fdr", n = 43)) %>% # correcting for 43 tests, not 68
        ungroup()
      partheritresults$annot.p <- if_else(partheritresults$fdr < 0.05, as.character(round(partheritresults$fdr, digits = 4)), "")
      partheritresults$significant = if_else(partheritresults$fdr < 0.05, "Yes", "")
      write.table(partheritresults, 
                  paste0("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/final_results/european_lr/results/results_tables/",unique(partheritresults$Annotation),"_results_FDR43.txt"),
                  sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
}

############
# Correct for 3 annotations
############

# SA

all_annots_results = read.table("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/final_results/european_lr/results_tables/all_annots.txt", 
                                fill = T, row.names = NULL)

all_annots_results = all_annots_results[-13]
colnames(all_annots_results) = c(colnames(all_annots_results)[2:12],"Analyis",colnames(all_annots_results)[13:16])

all_annots_results_left = all_annots_results[grep(all_annots_results$Region, pattern = "le_"),]
all_annots_results_right = all_annots_results[grep(all_annots_results$Region, pattern = "re_"),]
all_annots_results_left = all_annots_results_left[,-16]
all_annots_results_right = all_annots_results_right[,-16]

all_annots_results_left = all_annots_results_left %>% 
  group_by(Region) %>%
  mutate(fdr2 = my.p.adj(fdr, method = "fdr")) %>%
  ungroup()
all_annots_results_left$significant = if_else(all_annots_results_left$fdr2 < 0.05, "Yes", "")
write.xlsx(all_annots_results_left, 
           "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/final_results/european_lr/results_tables/all_annots_LE_FDR48_FDR3.xlsx",
           row.names = F, quote = F)

all_annots_results_right$fdr2 = all_annots_results_right %>% 
  group_by(Region) %>%
  mutate(fdr2 = p.adjust(fdr, method = "fdr")) %>%
all_annots_results_right$significant = if_else(all_annots_results_right$fdr2 < 0.05, "Yes", "")
write.xlsx(all_annots_results_right, 
           "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/final_results/european_lr/results_tables/all_annots_RI_FDR48_FDR3.xlsx",
           row.names = F, quote = F)

# DTI
# read all annotations and merge
neanRA_results_dti = read.table("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/final_results/dti/nean_RA_hg19-rinker_et_al.sorted_results_FDR25.txt", 
                                    fill = T, row.names = NULL, header = T)
neanRA_results_dti = neanRA_results_dti[-48,]
neanDep_results_dti = read.table("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/final_results/dti/neanDepRegions_hg19.sorted_results_FDR25.txt", 
                                fill = T, row.names = NULL, header = T)
fHGE_results_dti = read.table("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/final_results/dti/fetal_hge_hg19.merged.sorted_results_FDR25.txt", 
                                 fill = T, row.names = NULL, header = T)

all_annots_dti = rbind(neanRA_results_dti, neanDep_results_dti, fHGE_results_dti)

colnames(all_annots_dti) = c("to_delete", 
                                     colnames(all_annots_dti)[1:ncol(all_annots_dti)-1])
all_annots_dti = all_annots_dti[,-1]

all_annots_dti_left = all_annots_dti[grep(all_annots_dti$Region, pattern = "left"),]
all_annots_dti_right = all_annots_dti[grep(all_annots_dti$Region, pattern = "right"),]

all_annots_dti_left = all_annots_dti_left %>% 
  group_by(Region) %>%
  mutate(fdr2 = p.adjust(Enrichment_p, method = "fdr")) %>%
  ungroup()
all_annots_dti_left$significant = if_else(all_annots_dti_left$fdr2 < 0.05, "Yes", "")

write.xlsx(all_annots_dti_left, 
           "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/final_results/dti/all_annots_results_LE_FDR25_FDR3.xlsx",
           row.names = F, quote = F)

all_annots_dti_right = all_annots_dti_right %>% 
  group_by(Region) %>%
  mutate(fdr2 = p.adjust(fdr, method = "fdr")) %>%
  ungroup()
all_annots_dti_right$significant = if_else(all_annots_dti_right$fdr2 < 0.05, "Yes", "")

write.xlsx(all_annots_dti_right, 
           "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/final_results/dti/all_annots_results_RI_FDR25_FDR3.xlsx",
           row.names = F, quote = F)