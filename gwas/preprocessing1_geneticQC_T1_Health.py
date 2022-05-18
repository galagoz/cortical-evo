
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified Nov 2021

@author: B. Molz
"""
""" 
Create first filter including: Genetic QC, Availability of T1's, Disoder status
Main things that need to be motified are
- create the neuro healthy list first
- Import these exclusion codes for ICD10, 09 and non cancer illness --> all in one excel file
- if paths change : change accordingly in the scirpt
- imaging parameters can be modifed, currently checking if freesurfer was run with both T1 and T2 scans, T1 not unusable
- NOTE for the replication we exclude the first UKB release, here ukb21288, which should remove all particpants scanned prior to 2018
- you end up with both a final boolean mask that sorts data fields of interest and a final EID list that can be used to e.g. create a list for imaging processing
"""

#Import stuff
import pandas as pd
#-------------------------------------------
# REPLICATION SAMPLE
#-------------------------------------------
#Load list of IDs that survived genetic QC
gen_filter_wb = pd.read_csv ('/data/workspaces/lag/workspaces/lg-ukbiobank/derived_data/genetic_data/snp/subset_imagingT1_40k/v1_white_british_ancestry/SQC/imagingT1_ind_list_sqc_postrelatedness.txt',header = None, names=['eid'])
#Load ID's of current sample
eid = pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/current_release/ukb43760/ukb43760.csv', sep =",", usecols =[0])
mask_current_release = eid['eid'].isin(gen_filter_wb['eid'])
#Load list of IDs and T1s that were in the April 2018 r
t1_2018 = pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/old_releases/ukb21288/ukb21288.csv', sep =",", usecols =[1244])
eid_2018 =pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/old_releases/ukb21288/ukb21288.csv', sep =",", usecols =[0])
#  Check who had imaging data in 2018
t1_2018_filter = t1_2018.notna()
# get the ID's of the subjects with imaging data
eid_2018_filtered = eid_2018[t1_2018_filter['20252-2.0']]
# filter by these ID's --> NOTE: We DONT want these people so if their IDs are in the current ID's we want a FALSE
subsample_2ndrelease = ~eid['eid'].isin(eid_2018_filtered['eid'])
# concat the genetic white british filter with above 
gen_qc_subsample_2ndrelease = pd.concat([mask_current_release, subsample_2ndrelease], axis =1)
# just take indivduals which had a true in both e.g. they are not filter via genetic qc and are not part of the inital release
mask_replication = gen_qc_subsample_2ndrelease.all(axis =1)

#-------------------------------------------
# EUROPEAN SAMPLE
#-------------------------------------------
#Load list of IDs that survived genetic QC
gen_filter_we = pd.read_csv ('/data/workspaces/lag/workspaces/lg-ukbiobank/derived_data/genetic_data/snp/subset_imagingT1_40k/v2_white_ancestry/SQC/imagingT1_wa_ind_list_sqc_postrelatedness.txt',header = None, names=['eid'])
#Load ID's of current sample
mask_european = eid['eid'].isin(gen_filter_we['eid'])

#--------------------------------------------------------------------
#Load imaging related filter: 'Freesurfer run with T2 datafield'
#--------------------------------------------------------------------
mask_freesurfer =pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/current_release/ukb43760/ukb43760.csv', sep =",", usecols =['26500-2.0'])
#convert to boolean
mask_freesurfer = mask_freesurfer.fillna(0)
mask_freesurfer = mask_freesurfer.astype(bool)
#---------------------------------------------------------------------

#--------------------------------------------------------------------
# #Load health related filter: We only want neuro healthy particpants
#--------------------------------------------------------------------
mask_neurohealthy = pd.read_csv ('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/UKBB/ICD_diagnosis/ukb43760/mask_neuro_healthy_ukb43760.dat', header = None)
#---------------------------------------------------------------------
#also load the list of unusable T1's and fitler by these as the freesurfer output is likely corrupted
#--------------------------------------------------------------------
unusable_t1 = pd.read_fwf ('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/UKBB/unusable_T1/noT1_log.txt', header = None)
unusable_t1.columns =['eid']
mask_unusable = ~eid['eid'].isin(unusable_t1['eid'])
#---------------------------------------------------------------------

#-------------------------------------------
# REPLICATION SAMPLE
#-------------------------------------------
# now get overal maske and EID of current sample
enigma_evol_replication = pd.concat([mask_replication ,mask_freesurfer, mask_neurohealthy, mask_unusable], axis =1)
enigma_evol_replication = enigma_evol_replication.all(axis =1)
eid_enigma_evol_replication = eid[enigma_evol_replication]
#save 
enigma_evol_replication.to_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/replication/preprocessing/mask_enigmaEvol_replication_ukb43760.dat', index = False, header = False)
eid_enigma_evol_replication.to_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/replication/preprocessing/eid_enigmaEvol_replication_ukb43760.dat',header =False, index = False)


#-------------------------------------------
# European SAMPLE
#-------------------------------------------
# now get overal maske and EID of current sample
enigma_evol_european = pd.concat([mask_european ,mask_freesurfer, mask_neurohealthy, mask_unusable], axis =1)
enigma_evol_european = enigma_evol_european.all(axis =1)
eid_enigma_evol_european = eid[enigma_evol_european]
#save 
enigma_evol_european.to_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/preprocessing/mask_enigmaEvol_european_ukb43760.dat', index = False, header = False)
eid_enigma_evol_european.to_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/preprocessing/eid_enigmaEvol_european_ukb43760.dat',header =False, index = False)

