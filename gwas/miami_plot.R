#/bin/Rscript
#
# This script will generate a Miami plot.
#
##########
# Libraries
library(qqman)
library(data.table)
# Read data

gwas_replication = fread('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/replication/bgenieOutput/global/sumstats_ukb43760_global_surface_replication_allChr.gz', header = T, select= c(1,10,11,7))
gwas_e3 = fread('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/ENIGMA3_mixed_se_wo_Mean_Full_SurfArea_20190429_noGC.txt.gz', header = T, select= c(1,10,11,7))

gwas_asym = fread('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/SECA_local_version/sumstats/GCST90010427_buildGRCh37.tsv.gz', header = T)
gwas_surface_le = fread('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/SECA_local_version/sumstats/sumstats_ukb43760_global_surface_le_european_allChr.gz', header = T)

# remove chr X from asym sumstats, because it doesn't exist in the UKB sumstats
gwas_asym = gwas_asym[gwas_asym$chromosome != "X",]
gwas_asym$chromosome = as.integer(gwas_asym$chromosome)

# MPI orange = #ED6B06
# MPI green = #00786A

# Plot & Save

# replication vs. enigma3
png("/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/results/replication/replicationVSenigma3_Miamiplot_w1000_h500.png", type="cairo", width=1000, height=500)
par(mfrow=c(2,1))
par(mar=c(0,5,5,3))
manhattan(gwas_replication, ylim = c(0,35), chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("dodgerblue3", "grey80"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL)
par(mar=c(5,5,1,3))
manhattan(gwas_e3, ylim = c(35,0), chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("dodgerblue3", "grey80"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, xlab = "", xaxt = "n")
dev.off()

# asymmetry (Sha et al.) and left hemispheric surface area
png("/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/results/asymVSleftFull_Miamiplot_w1000_h500_v3.png", type="cairo", width=1000, height=500)
par(mfrow=c(2,1))
par(mar=c(0,5,5,3))
manhattan(gwas_asym, ylim = c(0,50), chr = "chromosome", bp = "base_pair_location", p = "p_value", snp = "variant_id", col = c("grey80", "dodgerblue3"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, main = "Asymmetry GWAS (Sha et al.)")
par(mar=c(5,5,1,3))
manhattan(gwas_surface_le, ylim = c(50,0), chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("grey80", "dodgerblue3"), chrlabs = NULL, suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = NULL, main = "Full left hemispheric surface area (UK Biobank)", xlab = "", xaxt = "n")
dev.off()
