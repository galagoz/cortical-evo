#!/bin/sh

#copy files after subsetting and variant QC from cluster to workspace 

for sample in european replication; do 
	#-----------------------------------------------------------------------------------------------
	#set paths
	#---------------------------------------------------------------------------------------------
	#path to output folder
	workspace_output=/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/${sample}/snpQC
	#clusterfs location
	cluster_output=/data/clusterfs/lag/users/barmol/enigma_evo/output/
	#---------------------------------------------------------------------------------------------------
	#set up folders and move files
	if [ ! -d ${workspace_output} ] && echo "Directory DOES NOT exists - creating it ."; then
		mkdir -p ${workspace_output}
	fi
	if [ ! -d ${workspace_output}subsetting ] && echo "Directory DOES NOT exists - creating it."; then
		mkdir -p ${workspace_output}/subsetting
	fi
	if [ ! -d ${workspace_output}/variantQC ] && echo "Directory DOES NOT exists - creating it."; then
		mkdir -p ${workspace_output}/variantQC
	fi
	#-----------------------------------------------------------------------------------------------------
	#let's move stuff around
	cp ${cluster_output}/subset_${sample}/* ${workspace_output}/subsetting
	cp -R ${cluster_output}/vqc_${sample}/* ${workspace_output}/variantQC
	echo 'done copying files for ${sample}' 
done
