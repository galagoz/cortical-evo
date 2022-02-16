## This script converts txt files to Rdata format

inputDir="/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface/withGlob/"
outputDir="/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface/withGlob/rdata/"

for (i in dir(inputDir, pattern="allChr.txt", all.files=F, full.names=F)) {
  tmp_dir=gsub(" ","", paste(inputDir,i), fixed=T)
  mergedGR=read.table(tmp_dir, header=T, fill=T, stringsAsFactors=F)
  RdataDir=paste(outputDir,gsub("txt", "Rdata", i))
  save(mergedGR,file=gsub(" ","", RdataDir))
}
