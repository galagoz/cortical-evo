# ========================================================
# Making uniform barplots to show LDSC partitioned heritability
# enrichment values. This would be run after making the 
# results tables. 
# 
# ========================================================

library(tidyverse)
library(scales)
library(BSDA)
library(here)

# Load the ordering of the brain regions, so your bar plots can be ordered front to back
regionordering = read.csv(here("scripts", "1000Genomes_Phase3_Analysis", "plotting", "freesurfer_orderandcolor.csv"))


filesGlobal = Sys.glob(here("/data/clusterfs/lag/users/gokala/enigma-evol/partherit/results_tables/*_results_FDR35.txt"))
#filesGlobal=c("/data/clusterfs/lag/users/gokala/enigma-evol/partherit/results_tables/HAR_results_FDR35.txt","/data/clusterfs/lag/users/gokala/enigma-evol/partherit/results_tables/Sweeps_results_FDR35.txt","/data/clusterfs/lag/users/gokala/enigma-evol/partherit/results_tables/NeanDepleted_results_FDR35.txt","/data/clusterfs/lag/users/gokala/enigma-evol/partherit/results_tables/NeanSNPs_3col_results_FDR35.txt")
for (i in 1:length(filesGlobal)){
  dataGlobal <- read.table(filesGlobal[i], header = TRUE, sep = "\t")
  dataGlobal$Region <- factor(dataGlobal$Region, levels = regionordering$Region)
  annot <- unique(dataGlobal$Annotation)
  pSA <- ggplot(data = dataGlobal[dataGlobal$Analysis == "Surface Area",], mapping = aes(Region, Enrichment))+
    geom_bar(stat="identity", position = "dodge", fill = "saddlebrown")+
    geom_linerange(position = position_dodge(width = 0.9), aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error))+
    geom_text(aes(x = Region, y = Enrichment+(Enrichment_std_error+2), label = annot.p), position = position_dodge(width = 0.9), size = 3,  fontface=3, hjust = 0, vjust = 0, angle = 45)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(x = "Region", 
         y = "Enrichment", 
         title = annot)
  pSA
  #ggsave(paste0("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/partherit/plots/enrichment/barplot_",annot,"_SA_MA6.svg"), width = 6, height = 3.25, units = "in")
  ggsave(paste0("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/partherit/plots/enrichment/barplot_",annot,"_SA_MA6.pdf"), width = 6, height = 3.25, units = "in")
  }

for (i in 1:length(filesGlobal)){
  dataGlobal = read.table(filesGlobal[i],header = TRUE, sep = "\t")
  dataGlobal$Region = factor(dataGlobal$Region,levels=regionordering$Region)
  annot = unique(dataGlobal$Annotation)
  pTH <- ggplot(data = dataGlobal[dataGlobal$Analysis == "Thickness",], mapping = aes(Region, Enrichment))+
    geom_bar(stat="identity", position = "dodge", fill = "saddlebrown")+
    scale_y_continuous(limits = c(-500, 500), oob = squish)+
    geom_linerange(position = position_dodge(width = 0.9), aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error))+
    geom_text(aes(x = Region, y = Enrichment-(Enrichment_std_error+2), label = annot.p), position = position_dodge(width = 0.9), size = 3,  fontface=3, hjust = 0, vjust = 0, angle = 45)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(x = "Region", 
         y = "Enrichment", 
         title = annot)
  #ggsave(paste0("./partherit/plots/barplot_",annot,"_TH_MA6.svg"), width =6, height = 3.25, units = "in")
  ggsave(paste0("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/partherit/plots/enrichment/barplot_",annot,"_TH_MA6.pdf"), width =6, height = 3.25, units = "in")
}
