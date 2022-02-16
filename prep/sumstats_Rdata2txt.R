## This script converts ancestry regressed summary statistics Rdata files to txt format,
## so that you can munge them and run LDSC partitioned heritability.

library(GenomicRanges)
inputDir="/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface_ancreg/"

for (i in list.files(inputDir, pattern=".Rdata", all.files=F, full.names=F)) {
  tmp_dir=gsub(" ","", paste(inputDir,i), fixed=T)
  print(paste0("Converting ", tmp_dir))
  load(tmp_dir)
  GWAS = as.data.frame(mcols(mergedGR))
  write.table(GWAS,gsub("Rdata","txt", tmp_dir, fixed=T),quote=F,row.names=F)
}
