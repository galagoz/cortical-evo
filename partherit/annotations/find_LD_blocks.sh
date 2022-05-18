####################################################################
##################### Haplotype Block Finder #######################
####################################################################

module load plink/1.9b6
module load bedtools/2.29.2

# Identify haplotype blocks using default max block length (200kb)
plink --bfile /data/clusterfs/lag/users/gokala/enigma-evol/partherit/beds/1KG_p3/plink/1KG_phase3_GRCh38_allchr --blocks no-pheno-req

# Extract haplotype blocks from chromosomes with Neanderthal Depleted Regions
awk '$1 == 1' plink.blocks.det > chr1.blocks.det
awk '$1 == 2' plink.blocks.det > chr2.blocks.det
awk '$1 == 7' plink.blocks.det > chr7.blocks.det
awk '$1 == 8' plink.blocks.det > chr8.blocks.det
awk '$1 == 18' plink.blocks.det > chr18.blocks.det

# Extract chr. #, haplotype start pos. and haplotype end pos. columns from .det files
awk '{print $1,$2,$3}' chr1.blocks.det > chr1.blocks.bed
awk '{print $1,$2,$3}' chr2.blocks.det > chr2.blocks.bed
awk '{print $1,$2,$3}' chr7.blocks.det > chr7.blocks.bed
awk '{print $1,$2,$3}' chr8.blocks.det > chr8.blocks.bed
awk '{print $1,$2,$3}' chr18.blocks.det > chr18.blocks.bed

# Add "chr" prefix to chr numbers in blocks.bed files
awk '$1="chr"$1' chr1.blocks.bed > tmp && mv tmp chr1.blocks.bed
awk '$1="chr"$1' chr2.blocks.bed > tmp && mv tmp chr2.blocks.bed
awk '$1="chr"$1' chr7.blocks.bed > tmp && mv tmp chr7.blocks.bed
awk '$1="chr"$1' chr8.blocks.bed > tmp && mv tmp chr8.blocks.bed
awk '$1="chr"$1' chr18.blocks.bed > tmp && mv tmp chr18.blocks.bed

# Sort 128 SNPs-list overlapping in introgressed and depleted annots
# This list is created in R using GenomicRanges package (reorganize_annotations.Rmd)
sort -k1,1 -k2,2n overlapping_nean_SNPs.bed > overlapping_nean_SNPs.sorted.bed

# Convert the white space delimited text file to tab delimited format
awk -v OFS='\t' '{$1=$1; print}' overlapping_nean_SNPs.sorted.bed > tmp && mv tmp overlapping_nean_SNPs.sorted.bed
awk -v OFS='\t' '{$1=$1; print}' chr1.blocks.bed > tmp && mv tmp chr1.blocks.bed
awk -v OFS='\t' '{$1=$1; print}' chr2.blocks.bed > tmp && mv tmp chr2.blocks.bed
awk -v OFS='\t' '{$1=$1; print}' chr7.blocks.bed > tmp && mv tmp chr7.blocks.bed
awk -v OFS='\t' '{$1=$1; print}' chr8.blocks.bed > tmp && mv tmp chr8.blocks.bed
awk -v OFS='\t' '{$1=$1; print}' chr18.blocks.bed > tmp && mv tmp chr18.blocks.bed

# Find intersects of Neanderthal Introgressed SNPs and 1KG haplotypes
intersectBed -wa -a chr1.blocks.bed \
		 -b overlapping_nean_SNPs.sorted.bed > chr1.introgression.blocks.bed
intersectBed -wa -a chr2.blocks.bed \
                 -b overlapping_nean_SNPs.sorted.bed > chr2.introgression.blocks.bed
intersectBed -wa -a chr7.blocks.bed \
                 -b overlapping_nean_SNPs.sorted.bed > chr7.introgression.blocks.bed
intersectBed -wa -a chr8.blocks.bed \
                 -b overlapping_nean_SNPs.sorted.bed > chr8.introgression.blocks.bed
intersectBed -wa -a chr18.blocks.bed \
                 -b overlapping_nean_SNPs.sorted.bed > chr18.introgression.blocks.bed

# Some SNPs in introgression-depletion overlap list fall within the same
# haplotype block. Thus, there are some duplicated rows in introgression.blocks.bed files. Extra rows are removed here.
sort chr1.introgression.blocks.bed | uniq > tmp && mv tmp chr1.introgression.blocks.bed
sort chr2.introgression.blocks.bed | uniq > tmp && mv tmp chr2.introgression.blocks.bed
sort chr7.introgression.blocks.bed | uniq > tmp && mv tmp chr7.introgression.blocks.bed
sort chr8.introgression.blocks.bed | uniq > tmp && mv tmp chr8.introgression.blocks.bed
sort chr18.introgression.blocks.bed | uniq > tmp && mv tmp chr18.introgression.blocks.bed

# Merge and sort introgressed haploypes from each chromosome
cat *.introgression.blocks.bed >> all.introgression.blocks.bed
sort -k1,1 -k2,2n all.introgression.blocks.bed > all.introgression.blocks.sorted.bed

# Subtract introgressed Neanderthal haplotypes from depleted regions
subtractBed -a neanDepRegions.sorted.bed -b all.introgression.blocks.sorted.bed > filtered.neanDepRegions.bed
