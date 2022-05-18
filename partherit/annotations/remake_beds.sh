#!/bin/bash
####################################################################
               # Merge Fetal Human Gained Enhancers #
# This script concatenates bed files, sorts rows based on chr. number
# and start pos., and merges overlapping ranges/rows.
####################################################################

module load bedtools/2.29.2

# Concatenation of bed files
cat 7pcw_Hu_gain_merged.bed 8_5pcw_Hu_gain_merged.bed 12Fpcw_Hu_gain_merged.bed 12Opcw_Hu_gain_merged.bed > fetal_hge.bed

# Sort concatenated bed file based on 1st and 2nd columns
sort -k1,1 -k2,2n fetal_hge.bed > fetal_hge.sorted.bed

# Remove the 4th column (placeholder) as it is obsolete
awk '{$4=""; print $0}' fetal_hge.sorted.bed > tmp && mv tmp fetal_hge.sorted.bed

# Convert the white space delimited text file to tab delimited format
awk -v OFS='\t' '{ $1=$1; print }' fetal_hge.sorted.bed > tmp && mv tmp fetal_hge.sorted.bed

# Check the total length of regions before merging overlapping regions
cat fetal_hge.sorted.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
wc -l fetal_hge.sorted.bed

# Merge overlapping rows/ranges
mergeBed -i fetal_hge.sorted.bed > fetal_hge.sorted.merged.bed

# Check the total length after merging overlapping regions
cat fetal_hge.sorted.merged.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
wc -l fetal_hge.sorted.merged.bed
chmod 777 fetal_hge.sorted.merged.bed # give read, write, execute permissions to any user

########################################################################
# Filter out Neanderthal introgressed SNPs from Nean. Depleted Regions #
# Second part of script kicks out Nean. introgressed SNPs out of 
# Nean. Depleted regions, which shouldn't be there anyway. But we
# found 128 such SNPs among ~123k Nean. introgressed SNPs.
########################################################################

# Remove the 4th column (placeholder) as it is obsolete
awk '{$4=""; print $0}' NeanDepleted_Vernot_Science_2016_TableS8.bed > neanDepRegions.bed

# Convert the white space delimited text file to tab delimited format
awk -v OFS='\t' '{ $1=$1; print }' neanDepRegions.bed > tmp && mv tmp neanDepRegions.sorted.bed

# Check the total length of regions before filtering out introgressed SNPs
cat neanDepRegions.sorted.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'

# Subtract introgressed SNPs from depleted regions
subtractBed -a neanDepRegions.sorted.bed -b NeanSNPs_3col.bed

# Check the total length of regions after filtering
cat neanDepRegions.sorted.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
