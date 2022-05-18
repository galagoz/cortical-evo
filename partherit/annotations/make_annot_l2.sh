#!/bin/bash
#########################################################
# Make .annot and .l2 files for LDSC part. heritability
# Download the latest version of bedtools to your local
# directory. Set path to that dir. with:
# export PATH="/home/gokala/bedtools2/bin:$PATH"
# Then you will not have any issues.
# You're welcome future me! :)
#########################################################

# Set variables
annot=$1
ldscDir="/home/gokala/ldsc/"
hapMap="/data/workspaces/lag/shared_spaces/Resource_DB/LDscores/Phase3/1000G_EUR_Phase3_baseline/print_snps.txt"
oneKg="/data/clusterfs/lag/users/gokala/enigma-evol/partherit/1000G_EUR_Phase3_plink/"
mainDir="/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/evolution/resources/new_annotations/beds/GRCh37/"
outDir="/data/clusterfs/lag/users/gokala/enigma-evol/partherit/new_annotations/nean_RA_hg19.sorted/"

#########################################################
# Make .annot files
# make_annot.py from LDSC is used
#########################################################

tmp_base=$(basename "$annot")
tmp_annot="${tmp_base%.*}"
#echo "##########################################"
echo "Making .annot files for $tmp_annot"
echo $annot
echo "##########################################"
mkdir ${outDir}${tmp_annot}

for chr in {1..22}; do
	echo $chr
	python ${ldscDir}make_annot.py \
		--bed-file $annot \
		--bimfile ${oneKg}1000G.EUR.QC.${chr}.bim \
		--annot-file ${outDir}${tmp_annot}/${tmp_annot}.${chr}.annot.gz
done
echo "$tmp_annot .annot files are ready!"
echo "##############"

#########################################################
# Make l2.ldscore files
# ldsc.py from LDSC is used
#########################################################

tmp_base=$(basename "$annot")
tmp_annot="${tmp_base%.*}"

echo "##########################################"
echo "Making .l2.ldscore files for $tmp_annot"
echo $annot
echo "##########################################"

for chr in {1..22}; do
	echo $chr
	python ${ldscDir}ldsc.py \
		--l2 \
		--bfile ${oneKg}1000G.EUR.QC.${chr} \
		--ld-wind-cm 1 \
		--print-snps $hapMap \
		--thin-annot \
		--annot ${outDir}${tmp_annot}/${tmp_annot}.${chr}.annot.gz \
		--out ${outDir}${tmp_annot}/${tmp_annot}.${chr}
	done

echo "$tmp_annot .l2.ldscore files are ready!"
echo "##############"
