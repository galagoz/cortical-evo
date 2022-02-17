## The folder containing GWAS summary statistics
fileloc = "/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface_ancreg/munged"
outputdir = "/data/clusterfs/lag/users/gokala/enigma-evol/final_analysis/data/replication/surface_ancreg/munged/"
filename = "munged_sumstats_list.txt"
curPattern = ".gz"

## read in gwas statistics file (compiled for all traits)
fGWASsumstats=gsub(" ","",paste(outputdir,filename))
write.table(dir(fileloc, pattern=curPattern, all.files=T,full.names=T),file=fGWASsumstats,row.names=F,col.names=F,quote=F)
