# ========================================================
# Making uniform barplots to show LDSC partitioned heritability
# enrichment values. This would be run after making the 
# results tables. 
# 
# ========================================================

library(tidyverse)
library(cowplot)
library(gridExtra)

# Load the ordering of the brain regions, so your bar plots can be ordered front to back
annots = list.files(path = "/data/clusterfs/lag/users/gokala/enigma-evol/dti/data/munged/results_tables/",
                    full.names = F, recursive = F, pattern = "FDR25")

plot_list = list()
#i=2
for (i in 1:length(annots)){

  # data wrangling  
  annotname = str_split(annots[i], pattern = "\\.")[[1]][1]
  resultsFile <- paste0("/data/clusterfs/lag/users/gokala/enigma-evol/dti/data/munged/results_tables/",annots[i])
  results <- read.table(resultsFile[1], header = TRUE, sep = "\t")
  
  results$Region = unlist(sapply(strsplit(results$Region, split = "_allChr_"), `[`, 1))
  
  results_left = results[grep(results$Region, pattern = "left"),]
  results_right = results[grep(results$Region, pattern = "right"),]
  results_global = results[-c(grep(results$Region, pattern = "left"), grep(results$Region, pattern = "right")),]
  
  # plot global meaures
  y_max1 <- max(results_global[results_global$Analysis == "WMtracts", 5])
  y_axis_max1 <- y_max1 + results_global[results_global$Enrichment == y_max1, 6] + 0.02
  label.df1 <- data.frame(Region = results_global$Region[results_global$significant=="Yes"&results_global$Analysis == "WMtracts"], 
                          Enrichment=results_global$Enrichment[results_global$significant=="Yes"&results_global$Analysis == "WMtracts"]+0.05)
  
  pWM_global = ggplot(data = results_global[results_global$Analysis == "WMtracts", ], mapping = aes(Region, Enrichment)) +
    geom_bar(stat = "identity", position = "dodge", fill = "#ED6B06") + ##00786A
    geom_linerange(position = position_dodge(width = 0.9), aes(ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "bottom") +
    labs(
      x = "Region",
      y = expression(paste("SNP-", italic("h"^{
        2
      }), " Enrichment")),
      title = "Global FA measures",
      subtitle = paste0("SNP-heritability enrichment in ",annotname)
    )
  
  pWM_global = pWM_global + geom_text(data = label.df1, label = "*",
                                      vjust = ifelse(label.df1$Enrichment < 0, 0.5, -0.5)) #,size=5
  
  if (any(results_global$Enrichment_p < 0.05) & any(!grepl(results_global$significant[which(results_global$Enrichment_p < 0.05)],pattern = "Yes"))) {
    label.df2 <- data.frame(Region = results_global$Region[results_global$Enrichment_p < 0.05], 
                            Enrichment = results_global$Enrichment[results_global$Enrichment_p < 0.05] + 0.05)
    pWM_global = pWM_global + geom_text(data = label.df2, label = "o",
                                        vjust = ifelse(label.df2$Enrichment < 0, 0.5, -0.5))
  }

  # plot left hemispheric meaures
  y_max1 <- max(results_left[results_left$Analysis == "WMtracts", 5])
  y_axis_max1 <- y_max1 + results_left[results_left$Enrichment == y_max1, 6] + 0.02
  label.df1 <- data.frame(Region = results_left$Region[results_left$significant=="Yes"&results_left$Analysis == "WMtracts"],
                          Enrichment=results_left$Enrichment[results_left$significant=="Yes"&results_left$Analysis == "WMtracts"]+0.05)
  
  pWM_left = ggplot(data = results_left[results_left$Analysis == "WMtracts", ], mapping = aes(Region, Enrichment)) +
    geom_bar(stat = "identity", position = "dodge", fill = "#ED6B06") + ##00786A
    geom_linerange(position = position_dodge(width = 0.9), aes(ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "bottom") +
    labs(
      x = "Region",
      y = expression(paste("SNP-", italic("h"^{
        2
      }), " Enrichment")),
      title = "Left hemispheric FA measures",
      subtitle = paste0("SNP-heritability enrichment in ",annotname)
    )
  
  pWM_left = pWM_left + geom_text(data = label.df1, label = "*",
                                  vjust = ifelse(label.df1$Enrichment < 0, 0.5, -0.5)) #,size=5
  if (any(results_left$Enrichment_p < 0.05) & any(!grepl(results_left$significant[which(results_left$Enrichment_p < 0.05)],pattern = "Yes"))) {
    label.df2 <- data.frame(Region = results_left$Region[as.numeric(results_left$Enrichment_p) < 0.05],
                            Enrichment=results_left$Enrichment[as.numeric(results_left$Enrichment_p) < 0.05]+0.05)
    pWM_left = pWM_left + geom_text(data = label.df2, label = "o",
                                    vjust = ifelse(label.df2$Enrichment < 0, 0.5, -0.5))
  }
  
  # plot right hemispheric meaures
  y_max1 <- max(results_right[results_right$Analysis == "WMtracts", 5])
  y_axis_max1 <- y_max1 + results_right[results_right$Enrichment == y_max1, 6] + 0.02
  label.df1 <- data.frame(Region = results_right$Region[results_right$significant=="Yes"&results_right$Analysis == "WMtracts"],
                          Enrichment=results_right$Enrichment[results_right$significant=="Yes"&results_right$Analysis == "WMtracts"]+0.05)
  
  pWM_right = ggplot(data = results_right[results_right$Analysis == "WMtracts", ], mapping = aes(Region, Enrichment)) +
    geom_bar(stat = "identity", position = "dodge", fill = "#ED6B06") + ##00786A
    geom_linerange(position = position_dodge(width = 0.9), aes(ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "bottom") +
    labs(
      x = "Region",
      y = expression(paste("SNP-", italic("h"^{
        2
      }), " Enrichment")),
      title = "Right hemispheric FA measures",
      subtitle = paste0("SNP-heritability enrichment in ",annotname)
    )
  
  pWM_right = pWM_right + geom_text(data = label.df1, label = "*",
                                    vjust = ifelse(label.df1$Enrichment < 0, 0.5, -0.5)) #,size=5
  if (any(results_right$Enrichment_p < 0.05) & any(!grepl(results_right$significant[which(results_right$Enrichment_p < 0.05)],pattern = "Yes"))) {
    label.df2 <- data.frame(Region = results_right$Region[as.numeric(results_right$Enrichment_p) < 0.05],
                            Enrichment=results_right$Enrichment[as.numeric(results_right$Enrichment_p) < 0.05]+0.05)
    pWM_right = pWM_right + geom_text(data = label.df2, label = "o",
                                      vjust = ifelse(label.df2$Enrichment < 0, 0.5, -0.5))
  }
  
  plot_list[[i]] = plot_grid(pWM_global, pWM_left, pWM_right,
                             labels = c("A", "B", "C"),
                             ncol = 1, nrow = 3)
}

ggsave(
  filename = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/partitioned_heritability/final_results/dti/WM_h2_enrichment_barplots_FDR25.pdf", 
  plot = marrangeGrob(plot_list, nrow=1, ncol=1), 
  width = 9, height = 15
)
