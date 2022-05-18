library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

fcorvals = "/data/clusterfs/lag/users/gokala/enigma-evol/sds/ancreg_replication/SDS_BJK_ancreg_replication_1kblocks_globals.csv"
fjason_corvals = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/sds/SDS_bjk_ancreg_1kblocks_fromJason.csv"
corvals = read.csv(fcorvals,stringsAsFactors=F)
jason_corvals = read.csv(fjason_corvals,stringsAsFactors=F)

#replication_thickness = corvals %>%
#  filter(grepl("(?=.*thickness)(?=.*globalCov)|(?=.*Thickness)(?=.*globalValues)",corvals$X,perl=TRUE)) %>%
#  select(X, global_corr_spearman)

#replace(replication_thickness$X,1,"Mean_Full_Thickness")
#replication_thickness$mergeCol=sapply(replication_thickness$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})

replication_surface = corvals %>%
  filter(grepl("(?=.*surface)(?=.*withGlob)",corvals$X,perl=TRUE)) %>%
  #filter(grepl("(?=.*surfaceDK)(?=.*global)|(?=.*averaged)(?=.*surface)",corvals$X,perl=TRUE)) %>%
  select(X, global_corr_spearman, BJK_ESTIM_Z)

#replication_surface$X[1]="Mean_Full_SurfArea"
replication_surface$mergeCol=sapply(replication_surface$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})
names(replication_surface)[names(replication_surface) == "X"] <- "rep_surf"
#replication_surface$color="red"

jason_corvals_surface = jason_corvals %>%
  filter(grepl("(?=.*surfavg)(?=.*Mean)|(?=.*SurfArea)(?=.*Full)",jason_corvals$X,perl=TRUE)) %>%
  select(X, global_corr_spearman, BJK_ESTIM_Z)

jason_corvals_surface$mergeCol=sapply(jason_corvals_surface$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})

merged = merge(jason_corvals_surface,replication_surface,by="mergeCol")

##########################################################
# MS version
# Correlation plot
ggscatter(merged, x = "global_corr_spearman.x", y = "global_corr_spearman.y", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Tilot et al. (2021)", ylab = "UK Biobank Replication", title = "tSDS-GWAS effect size correlations per cortical region",
          color = "#EEA300")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/plots/sds/e3_vs_replication_SDS_corcoef_scatter.pdf",
       width = 8, height = 8)

# Z-score plot
ggscatter(merged, x = "BJK_ESTIM_Z.x", y = "BJK_ESTIM_Z.y", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Tilot et al. (2021)", ylab = "UK Biobank Replication", title = "tSDS-GWAS effect size Z-scores per cortical region",
          color = "#49A7DE")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/results/plots/sds/e3_vs_replication_SDS_zscore_scatter_blue.pdf",
       width = 8, height = 8)

#########################################################

#pdf("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/surface_e3_vs_replication_scatter2.pdf")
#plot(merged$global_corr_spearman.x,merged$global_corr_spearman.y,xlab="ENIGMA E3",ylab="Replication")
#dev.off()

mergedTall = merged %>% gather(key=Data,value=CorCoef,global_corr_spearman.x,global_corr_spearman.y)
#ggplot(mergedTall, aes(mergeCol, CorCoef, fill = Data)) + geom_col(position = "dodge") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(mergedTall, aes(mergeCol, CorCoef, fill = Data)) + geom_col(position = "dodge") + coord_flip() + scale_fill_discrete(name = "Source", labels = c("E3", "Replication")) + ggtitle("SDS correlation coefficient comparison (surface)") + ylab("Correlation Coefficient") + xlab("Region")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/surface_e3_vs_replication.pdf",height = 10,width = 7)

#### Read-in left-right SDS corcoef.s

fcorvals_lr = "/data/clusterfs/lag/users/gokala/enigma-evol/sds/SDS_bjk_ancreg_1kblocks_regional_with_global.csv"
corvals_lr = read.csv(fcorvals_lr,stringsAsFactors=F)

left = corvals_lr %>%
  filter(grepl("(?=.*le_)",corvals_lr$X,perl=TRUE)) %>%
  #filter(grepl("(?=.*_le)(?=.*surfaceDK)",corvals_lr$X,perl=TRUE)) %>%
  select(X, global_corr_spearman)
left$mergeCol=sapply(left$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})
merged_left=merge(replication_surface,left,by="mergeCol")
#left$X[1]="Mean_Full_SurfArea"
names(left)[names(left) == "X"] <- "left"
left$color="green"

right = corvals_lr %>%
  filter(grepl("(?=.*re_)",corvals_lr$X,perl=TRUE)) %>%
  #filter(grepl("(?=.*_re)(?=.*surfaceDK)",corvals_lr$X,perl=TRUE)) %>%
  select(X, global_corr_spearman)
right$mergeCol=sapply(right$X,function (x) {unlist(strsplit(x,"_",fixed=TRUE))[2]})
merged_right=merge(replication_surface,right,by="mergeCol")
#right$X[1]="Mean_Full_SurfArea"
names(right)[names(right) == "X"] <- "right"
right$color="blue"

mergedTallLeft = merged_left %>% gather(key=Data,value=CorCoef,global_corr_spearman.x,global_corr_spearman.y)
ggplot(mergedTallLeft, aes(mergeCol, CorCoef, fill = Data))+ coord_flip()  + geom_col(position = "dodge") + scale_fill_discrete(name = "Source", labels = c("Left-Right Averaged", "Left Hemisphere only")) + ggtitle("SDS correlation coefficient comparison (surface areas)") + ylab("Correlation Coefficient") + xlab("Region")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/surfarea_leftRightAveraged_vs_left.pdf",height = 10,width = 7)

mergedTallRight = merged_right %>% gather(key=Data,value=CorCoef,global_corr_spearman.x,global_corr_spearman.y)
ggplot(mergedTallRight, aes(mergeCol, CorCoef, fill = Data)) + geom_col(position = "dodge") + coord_flip() + scale_fill_discrete(name = "Source", labels = c("Left-Right Averaged", "Right Hemisphere only")) + ggtitle("SDS correlation coefficient comparison (surface areas)") + ylab("Correlation Coefficient") + xlab("Region")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/surfarea_leftRightAveraged_vs_right.pdf",height = 10,width = 7)

# Final plot

df_list=list(replication_surface,left,right)
all_merged=Reduce(function(x, y) merge(x, y, all=TRUE, by="mergeCol"), df_list, accumulate=F)
#all_merged = merge(jason_corvals_surface,replication_surface,left,right,by="mergeCol")

mergedTall = all_merged %>% gather(key=Data,value=CorCoef,global_corr_spearman.x,global_corr_spearman.y,global_corr_spearman)
ggplot(mergedTall, aes(reorder(mergeCol,CorCoef), CorCoef, fill = Data)) + geom_col(position = "dodge") + coord_flip() + scale_fill_discrete(name = "Data", labels = c("Left Hemisphere", "European subset", "Right Hemisphere")) + ggtitle("SDS correlation coefficients (Surface Area)") + ylab("Correlation Coefficient") + xlab("Region")
#ggplot(mergedTall, aes(mergeCol, CorCoef, fill = Data)) + geom_col(position = "dodge") + coord_flip() + scale_fill_discrete(name = "Source", labels = c("E3", "Replication", "Left Hemisphere", "Right Hemisphere")) + ggtitle("SDS correlation coefficients (Surface Area)") + ylab("Correlation Coefficient") + xlab("Region")
ggsave("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evol-pipeline/eurSubset_and_l-r_w-globalCov_surfArea.pdf",height = 10,width = 7)