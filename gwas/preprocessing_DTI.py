#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 16:46:50 2021

@author: B Molz
"""

#This script prepares phenotypes and covariates for associtatoin analysis

"""Several filtering and QC steps are applied, filter masks and overview plots created. As last steps files are reformatted to BGENIE standards 
"""

#Import stuff
import pandas as pd
import os
import numpy as np
from statsmodels import robust
import matplotlib.pyplot as plt
import seaborn as sns
import csv

#-------------------------------------------
# create directory for DTI output
#------------------------------------------
pathDTI='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/DTI' 
if not os.path.exists(pathDTI):
    os.mkdir(pathDTI)
#------------------------------------------------------------------------------
# 1 Get the overall sample
#------------------------------------------------------------------------------

#Load list of IDs that survived genetic QC
gen_filter_wb = pd.read_csv ('/data/workspaces/lag/workspaces/lg-ukbiobank/derived_data/genetic_data/snp/subset_imagingT1_40k/v1_white_british_ancestry/SQC/imagingT1_ind_list_sqc_postrelatedness.txt',header = None, names=['eid'])
#Load ID's of current sample
eid = pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/current_release/ukb43760/ukb43760.csv', sep =",", usecols =[0])
mask_current_release = eid['eid'].isin(gen_filter_wb['eid'])
# EUROPEAN SAMPLE
#Load list of IDs that survived genetic QC
gen_filter_we = pd.read_csv ('/data/workspaces/lag/workspaces/lg-ukbiobank/derived_data/genetic_data/snp/subset_imagingT1_40k/v2_white_ancestry/SQC/imagingT1_wa_ind_list_sqc_postrelatedness.txt',header = None, names=['eid'])
#Load ID's of current sample
mask_european = eid['eid'].isin(gen_filter_we['eid'])
# #Load health related filter: We only want neuro healthy particpants
mask_neurohealthy = pd.read_csv ('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/UKBB/ICD_diagnosis/ukb43760/mask_neuro_healthy_ukb43760.dat', header = None)
#also load the list of unusable T1's and fitler by these as the freesurfer output is likely corrupted
unusable_t1 = pd.read_fwf ('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/UKBB/unusable_T1/noT1_log.txt', header = None)
unusable_t1.columns =['eid']
mask_unusable = ~eid['eid'].isin(unusable_t1['eid'])
#---------------------------------------------------------------------
# now get overal maske and EID of current sample
enigma_evol_european_DTI = pd.concat([mask_european, mask_neurohealthy, mask_unusable], axis =1)
enigma_evol_european_DTI = enigma_evol_european_DTI.all(axis =1)
eid_enigma_evol_european_DTI = eid[enigma_evol_european_DTI]
#save 
enigma_evol_european_DTI.to_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/preprocessing/mask_enigmaEvol_europeanDTI_ukb43760.dat', index = False, header = False)
eid_enigma_evol_european_DTI.to_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/preprocessing/eid_enigmaEvol_europeanDTI_ukb43760.dat',header =False, index = False)

#---------------------------------------------------------
#2 Now we extract DTI TBSS mean FA values 
#--------------------------------------------------------------
DTI_TBSS =pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/current_release/ukb43760/ukb43760.csv',usecols=[0,*range(5290,5385,2)], dtype = object)
#now we filter the phenotype by availability of DTI metrics and previous exclusioin criteria
DTI_european = DTI_TBSS[enigma_evol_european_DTI].copy()
DTI_european_noNaN=DTI_european.dropna(how = 'any').copy()
cols=list(DTI_european_noNaN.iloc[:,1:])
DTI_european_noNaN[cols]=DTI_european_noNaN[cols].astype(float)
DTI_european_noNaN['eid']=DTI_european_noNaN.eid.astype(str)
#save as a CSV excel file 
DTI_european_noNaN.to_excel('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/DTI/european_DTITBSS.xlsx', index= None)

#--------------------------------------------------------------------------------
# 3 PHENOTYPE QUALITY CONTROL
#-----------------------------------------------------------------------------------
global_DTI=DTI_european_noNaN.iloc[:,1:7]
right_DTI=DTI_european_noNaN.iloc[:,7:48:2]
left_DTI=DTI_european_noNaN.iloc[:,8:49:2]
lower_bound =[]
upper_bound =[]
for column in global_DTI.iloc[:,0:]:
    lower = global_DTI[column].median() - 5*robust.mad(global_DTI[column])
    upper= global_DTI[column].median() + 5*robust.mad(global_DTI[column])
    lower_bound.append(lower)
    upper_bound.append(upper)
##empty dataframe
global_DTI_qc = pd.DataFrame()   
#check each column for outliers, make qc'ed dataframe
for ind, column in enumerate(global_DTI.iloc[:,0:].columns):
    global_DTI_qc[column] = global_DTI[column].mask((global_DTI[column] < lower_bound[ind]) | (global_DTI[column] > upper_bound[ind]))

# do the same thing for left right but here we keep the N similar across hemispheres
lofdfs = [right_DTI, left_DTI]
rep_pheno_qc = {}

for key,df in enumerate(lofdfs):
#QC each by using the IQR 
    lower_bound =[]
    upper_bound =[]
    for column in df.iloc[:,0:]:
        lower =df[column].median() - 5*robust.mad(df[column])
        upper= df[column].median() + 5*robust.mad(df[column])
        lower_bound.append(lower)
        upper_bound.append(upper)
        qc=pd.DataFrame()   
    for ind, column in enumerate(df.iloc[:,0:].columns):
        qc[column] = df[column].mask((df[column] < lower_bound[ind]) | (df[column] > upper_bound[ind]))
        #also immediatelly delete the global outliers
        rep_pheno_qc[key]=qc

right_DTI_qc=rep_pheno_qc[0]
left_DTI_qc=rep_pheno_qc[1]

# quirk of pandas - can't work on columns in a global way otherwise so we save them to a list 
dtiLE_col =left_DTI_qc.columns.values.tolist()
dtiRE_col =right_DTI_qc.columns.values.tolist()
#and replace them with numbers
left_DTI_qc.columns=np.arange(len(left_DTI_qc.columns))
right_DTI_qc.columns=np.arange(len(right_DTI_qc.columns))

#now we replace outliers with NAN by checking each dataframe frame and get sumstats
dtiLE_masked = left_DTI_qc.mask(right_DTI_qc.isna()).copy()
dtiRE_masked = right_DTI_qc.mask(left_DTI_qc.isna()).copy()

#now name the columns again
dtiLE_masked.columns = dtiLE_col
dtiRE_masked.columns = dtiRE_col

#get all the QCed data in one dataframe, add eid for later use
final_DTI_qc =pd.concat([DTI_european_noNaN['eid'],global_DTI_qc,dtiLE_masked,dtiRE_masked],axis=1)
dtifinalsumstats = final_DTI_qc.iloc[:,1:].describe() 
dtifinalsumstats.loc['percentange',:] =dtifinalsumstats.loc['count',:] / len(DTI_european_noNaN.index)
dtifinalsumstats.to_csv('{}/ukb43760_dtiTBSS_QCed_sumstats.csv'.format(pathDTI))

#save EIDs for later use
final_DTI_qc['eid'].to_csv('{}/eid_ukb43760_dtiTBSS.dat'.format(pathDTI),index=None,header=None)
#get the overall ID list to create final mask
eid = pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/current_release/ukb43760/ukb43760.csv', sep =",", usecols =[0])
mask_QCed = eid['eid'].isin(final_DTI_qc['eid'])
#let's print out the N just as a sanity check
mask_QCed.sum()
#save the final mask - this combines neuro healthy, imaging, pheno QC and genetic sample QC
mask_QCed.to_csv('{}/mask_ukb43760_dtiTBSS.dat'.format(pathDTI),index=None,header=None)

final_DTI_qc.to_excel('{}/ukb43760_dtiTBSS_qced.xlsx'.format(pathDTI),sheet_name='ukb43760_dtiTBSS_FA')

# quick check if all used EIDs have been used for surface based metrics --> do we need to rerun subsetting?
eid_surfacebased=pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/preprocessing/eid_afterQC_european_ukb43760.dat',header=None,dtype=str)
eid_surfacebased.columns=['eid']
#chekch the overlap 
overlap=final_DTI_qc['eid'].isin(eid_surfacebased['eid'])
print ('the differnce between this and the original sample is',len(final_DTI_qc) - overlap.sum())

#---------------------------------------------------------------------------
# 4 PLOT all regional metrics after QC as well 
#--------------------------------------------------------------------------
# make violinplot with Seaborn
# loop over it to make it more readable 
plt.figure(figsize=(30,7))
bplot=sns.violinplot(y='value', x='variable', 
                 data=pd.melt(final_DTI_qc.iloc[:,1:]), 
                 width=0.5,
                 palette="GnBu_d")
bplot.set_xticklabels(bplot.get_xticklabels(), rotation=45, horizontalalignment='right')
sns.despine()
bplot.set_xlabel("ROI",fontsize=12)
bplot.set_ylabel("mean FA",fontsize=12)
bplot.tick_params(labelsize=12)
bplot.autoscale()
fig = bplot.get_figure()
fig.savefig('{}/ukb43760_DTI_violinPlots.png'.format(pathDTI),bbox_inches="tight")

#---------------------------------------------------------------------------
# 5 Prep Covariates
#--------------------------------------------------------------------------
# get the other covariates used in the analysis 
confounds =pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/current_release/ukb43760/ukb43760.csv', sep =",", usecols =['eid','31-0.0','21003-2.0','22009-0.1','22009-0.2','22009-0.3','22009-0.4','22009-0.5','22009-0.6','22009-0.7','22009-0.8', '22009-0.9', '22009-0.10','54-2.0','22000-0.0'])
# sort the columns again as pandas rearranges them... 
confounds = confounds[['eid','31-0.0','21003-2.0','22009-0.1','22009-0.2','22009-0.3','22009-0.4','22009-0.5','22009-0.6','22009-0.7','22009-0.8', '22009-0.9', '22009-0.10','54-2.0','22000-0.0']]
confounds.columns = ['eid','sex','age_at_scan','PC1','PC2','PC3','PC4','PC5', 'PC6','PC7', 'PC8','PC09','PC10','assessmentcentre','genotypearray']
#mask the dataframe to fit current sample
confound_dti = confounds[mask_QCed].copy()
#change sex value from 0 female 1 male to female = 2 male =1 
confound_dti.loc[confound_dti['sex'] == 0, 'sex'] = 2
#make genotype batch into a binary array variables
confound_dti.loc[confound_dti['genotypearray'] < 0, 'genotypearray'] = 0
confound_dti.loc[confound_dti['genotypearray'] > 0, 'genotypearray'] = 1
# get mean age in sample
age_mean =confound_dti['age_at_scan'].mean()
# get mean centred age
confound_dti.loc[:,'ageCSq'] = confound_dti.loc[:,'age_at_scan'].apply(lambda x: np.square(x - age_mean))
# get sex age interaction
confound_dti.loc[:,'age_sex'] = confound_dti.loc[:,'age_at_scan'] * confound_dti.loc[:,'sex']
#get sex mean centred age interaction
confound_dti.loc[:,'ageCSq_sex'] =confound_dti.loc[:,'ageCSq'] * confound_dti.loc[:,'sex']
#let's factorize the categorical variables 
cat_variables=['sex','assessmentcentre','genotypearray']
for col in cat_variables:
    confound_dti[col] = pd.factorize(confound_dti[col])[0]
#also as we dummy code to this after subsetting the inital array as we otherwise end up with more dummy variables than needed
confound_dti=pd.get_dummies(confound_dti, columns=['assessmentcentre'],drop_first =True)
#save it
confound_dti.to_excel('{}/ukb43760_DTI_covariates.xlsx'.format(pathDTI),index = None)

#---------------------------------------------------------------------------
# 6 Prep for BGENIE
#--------------------------------------------------------------------------

pathDTI_gwas='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/DTI/GWAS' 
if not os.path.exists(pathDTI_gwas):
    os.mkdir(pathDTI_gwas)
    
#-----------------------
#read the sample file, strip first row ( 0 0 0) and only import one column
bgen_list = pd.read_csv('{}/snpStats/subsetting/ukb43760_enigmaEvo_europeanDTI_chr10.sample'.format(pathDTI),delim_whitespace=True, usecols=[0],skiprows=[1],dtype=str)
# read the roi list - needed later to give better column names (currently UKB IDPs)
#roi_name = pd.read_excel('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/QC_5MAD/ROI_table.xlsx', header =None)
#--------------------------------
# phenotypes
#--------------------------------
DTI_gwas=final_DTI_qc.set_index('eid')
# sort by the sample file
DTI_gwas_sorted = DTI_gwas.loc[bgen_list['ID_1']]
DTI_gwas_sorted.fillna(-999, inplace=True)
DTI_gwas_sorted.to_csv('{}/ukb43670_DTI_phenotype_bgenie.table'.format(pathDTI_gwas), index = None, sep=' ', quoting = csv.QUOTE_NONE)
    
#--------------------------
# covariates
#--------------------------
confound_dti['eid']=confound_dti.eid.astype(str)
confound_DTI_gwas=confound_dti.set_index('eid')
# sort by the sample file
confound_DTI_gwas_sorted = confound_DTI_gwas.loc[bgen_list['ID_1']]
confound_DTI_gwas_sorted.fillna(-999, inplace=True)
confound_DTI_gwas_sorted= confound_DTI_gwas.round(decimals=5)
confound_DTI_gwas_sorted.to_csv('{}/ukb43670_DTI_confound_bgenie.table'.format(pathDTI_gwas), index = None, sep=' ', quoting = csv.QUOTE_NONE)

#--------------------------
# SNPs to keep - original files created by variantQC procedure
#--------------------------
file_name = '{}/snpStats/variantQC/ukb43760_enigmaEvo_europeanDTI_chr{}.snpstats_mfi_hrc.snps2keep'
output_name = '{}/ukb43760_enigmaEvo_europeanDTI_snps2keep_chr{}.table'

for i in range(1, 23):
    current_data=pd.read_csv(file_name.format(pathDTI,i),sep='\s+',usecols =['RS_ID.UKB'])
    current_data.to_csv(output_name.format(pathDTI_gwas,i), index = None, sep=' ', quoting = csv.QUOTE_NONE)