# This script reads the raw bed files and
# saves each column as a field for further
# reformatting in bash.

options(stringsAsFactors=FALSE)

fbeds = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/raw_beds/"
beds = list.files(fbeds,pattern = "\\.")

for (i in beds) {
  tmp_bed = read.table(paste0(fbeds,i))
  tmp_bed = read.table(fbeds)
  #tmp_bed = tmp_bed[,c(1,2,3)]
  #colnames(tmp_bed) = c("chr", "start", "end")
  #tmp_bed$length = (tmp_bed$end - tmp_bed$start)
  #annot_lengths$length[which(annot_lengths$annots==i)] = sum(tmp_bed$length)
  write.table(tmp_bed,paste0("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/",i),quote=F,row.names=F,col.names=F,sep="\t")
}

# Here the annotation lenghts are measured.    # note from future-self: lenghts are not important, SNP prop.s are.
# Lenghts are important to see if annotations
# are large enough. However, the correct check
# is to see how many 1000G phase 3 genomes fall
# within your annotation (should be at least
# 1% of the whole genome).
fbeds = "/data/workspaces/lag/workspaces/lg-genlang/Working/Evolution/results/ldsc_cts_annots/cell_type_specific/"
beds = list.files("/data/workspaces/lag/workspaces/lg-genlang/Working/Evolution/results/ldsc_cts_annots/cell_type_specific/",pattern = "\\.")

annot_lengths=data.frame(annots=beds,length=0)

for (i in beds) {
  tmp_bed = read.table(paste0(fbeds,i))
  tmp_bed = tmp_bed[,c(1,2,3)]
  colnames(tmp_bed) = c("chr", "start", "end")
  tmp_bed$length = (tmp_bed$end - tmp_bed$start)
  annot_lengths$length[which(annot_lengths$annots==i)] = sum(tmp_bed$length)
}

write.table(annot_lengths,"/data/workspaces/lag/workspaces/lg-genlang/Working/Evolution/results/ldsc_cts_annots/annot_sizes.txt",quote=F,row.names=F)

for (i in beds) {
  tmp_bed = read.table(paste0(fbeds,i))
  #tmp_bed = tmp_bed[,c(1,2,3)]
  #colnames(tmp_bed) = c("chr", "start", "end")
  #tmp_bed$length = (tmp_bed$end - tmp_bed$start)
  #annot_lengths$length[which(annot_lengths$annots==i)] = sum(tmp_bed$length)
  write.table(tmp_bed,paste0("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/",i),quote=F,row.names=F,col.names=F,sep="\t")
}

# Quick fix for .bed files with non-tab 
# delimiter. Reads .beds and writes back
# to the same path with sep="\t".

bed1=read.table("/data/clusterfs/lag/users/gokala/enigma-evol/partherit/beds/new_beds2/neanderthal_NDA_and_RA_hg19-rinker_et_al.sorted.bed")
bed2=read.table("/data/clusterfs/lag/users/gokala/enigma-evol/partherit/beds/new_beds2/nean_introg_SNPs_hg19.sorted.bed")
write.table(bed1,"/data/clusterfs/lag/users/gokala/enigma-evol/partherit/beds/new_beds2/neanderthal_NDA_and_RA_hg19-rinker_et_al.sorted2.bed",row.names = F,
                                                                                                                                             col.names = F,
                                                                                                                                             quote = F,
                                                                                                                                             sep = "\t")
write.table(bed2,"/data/clusterfs/lag/users/gokala/enigma-evol/partherit/beds/new_beds2/nean_introg_SNPs_hg19.sorted2.bed",row.names = F,
                                                                                                                          col.names = F,
                                                                                                                          quote = F,
                                                                                                                          sep = "\t")