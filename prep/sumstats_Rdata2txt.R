## This script converts Rdata files to txt format
library(GenomicRanges)
inputDir="/data/clusterfs/lag/users/gokala/ancestry_regression/AncestryRegressionData_noGC/"
#"P:/workspaces/lg-genlang/Working/Evolution/AncestryRegressionData_noGC/"

for (i in list.files(inputDir, pattern="ancreg.Rdata", all.files=F, full.names=F)) {
  tmp_dir=gsub(" ","", paste(inputDir,i), fixed=T)
  print(paste0("Converting ", tmp_dir))
  load(tmp_dir)
  GWAS = as.data.frame(mcols(mergedGR))
  write.table(GWAS,gsub("Rdata","txt", tmp_dir, fixed=T))
}