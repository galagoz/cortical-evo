#!/bin/bash

# Gokberk Alagoz
# 29.04.2021
################

# This script extracts 1000G phase 3 SNPs which 
# fall within the evolutionary annotations from 
# VCF files.

################
# paths:
vcf="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/1000G_phase3"
bed="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/for_SNP_counts"
outDir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/1KGP3_snps"
plinkDir="/data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh38/plink"
plinkDirEUR="/data/workspaces/lag/shared_spaces/Resource_DB/1KG_phase3/GRCh37/plink"
plinkDirTilot="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/1KGP3_EUR_snps_fromAmanda"
sites="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz"
plinkDirEURnonFIN_GRCh38="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/1KGP3_EURnonFIN_subset"

# modules
module load bedtools/2.29.2
module load vcftools/0.1.17
module load plink/1.9b6

# 1) Sort your bed files based on chr and pos
# and then split by chromosome.

#for i in $bed/*.bed; do
#	echo "Sorting "${i}
#	sort -k1,1 -k2,2n ${i} > ${i%.bed}_sorted.bed;
#done

#for i in $bed/*sorted.bed; do
#	for j in `bedextract --list-chr ${i}`; do
#		bedextract ${j} ${i} > ${i%.bed}.${j}.bed;
#       	done;
#done

# 2) Filter VCF files using bed files. Then you
# can see how many 1KGP3 SNPs your annotation
# harbours.

# 2.a) Used vcftools here, but the code is taking ages
# to run!!!

#for i in $bed/*sorted.bed; do
#	tmp_bed=$(basename "$i")
#
#	echo "###############################"
#	echo "Working on "$tmp_bed
#	echo "###############################"
#
#	for j in {1..22}; do
#	echo "Filtering chromosome "$j
#	vcftools --gzvcf ${vcf}/ALL.chr${j}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
#	       	 --bed ${i%.bed}.chr${j}.bed --recode \
#		 --stdout | gzip -c > ${outDir}/${tmp_bed%.bed}_chr${j}_1KGP3_filtered.vcf.gz;
#	done;
#done

# 2.b) Converted Amanda's python script to bash:
# (Step0_find_SNPs_in_1000G_Phase3.py)
# Hopefully it is faster!

# Problem: "--extract range" flag requires bed
# files with 4 columns. Last column being "arbitrary
# range IDs". I have no clue what that means, so I 
# added a dummy code as 4th column (all 1s).
# for i in *.sorted.bed; do awk '{print $0,"\t"1}' $i > ${i%.bed}_forSNPCounts.bed; done
# *sorted_forSNPCounts.bed

#for i in $bed/*weiss_et_al_forSNPCounts.bed; do
#        tmp_bed=$(basename "$i")
#
#	echo "###############################"
#        echo "Working on "$tmp_bed
#        echo "###############################"
#	
#	plink --bfile ${plinkDir}/1KG_phase3_GRCh38_allchr \
#		--make-just-bim --allow-no-sex \
#		--extract range "${i}" --out ${outDir}/1KG_phase3_GRCh38_allchr.${tmp_bed}
#
#	echo "Done with "${i}
#	echo "###############################";
#done

# 2.c) Use EUR_nonFIN and see if a similar
# number of SNPs survive. It is GRCh37/hg19
# though - not sure if it makes sense.

#for i in $bed/*weiss_et_al_forSNPCounts.bed; do
#        tmp_bed=$(basename "$i")
#
#        echo "###############################"
#        echo "Working on "$tmp_bed
#        echo "###############################"
#
#        plink --bfile ${plinkDirEUR}/1KG_phase3_GRCh37_EUR_nonFIN_allchr \
#                --make-just-bim --allow-no-sex \
#                --extract range "${i}" --out ${outDir}/1KG_phase3_GRCh37_EUR_nonFIN_allchr.${tmp_bed}
#
#        echo "Done with "${i}
#        echo "###############################";
#done

# 2.d) Use EUR from Amanda and see if you
# have similar SNP numbers.

#for i in $bed/*.bed; do
#        tmp_bed=$(basename "$i")
#
#        echo "###############################"
#        echo "Working on "$tmp_bed
#        echo "###############################"
#
#        plink --bfile ${plinkDirTilot}/1000G.EUR.QC \
#		--make-just-bim --allow-no-sex \
#                --extract range "${i}" --out ${outDir}/1000G.EUR.QC_allchr.${tmp_bed}
#
#        echo "Done with "${i}
#        echo "###############################";
#done

# 2.e) I found another way to count SNPs
# in bed files. Try it and see if results
# match. Reference:
# https://www.biostars.org/p/9467869/#9468106

#for i in $bed/*.bed; do
#  echo -n "$i "
#  /home/gokala/bin/bcftools view -HR $i $sites | wc -l
#done

# 2.f) I made the EURnonFIN GRCh38 .bed subset myself.
# Now will run it with that one.

for i in $bed/*_forSNPCounts.bed; do
        tmp_bed=$(basename "$i")
	
        echo "###############################"
        echo "Working on "$tmp_bed
        echo "###############################"

        plink --bfile ${plinkDirEURnonFIN_GRCh38}/1KG_phase3_GRCh38_EUR_allchr \
                --make-just-bim --allow-no-sex \
                --extract range "${i}" --out ${outDir}/1KG_phase3_GRCh38_EURnonFIN_allchr.${tmp_bed%.bed}

        echo "Done with "${i}
        echo "###############################";
done
