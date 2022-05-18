#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified Nov 2022

@author: B. Molz
"""
'''
This script imports the created excel phenotype files and does some general QC:
Threshold +/- 5* MAD
If global outlier --> delete whole subject
If regionl outlier --> delete all metrics + both hemispheres e.g. same N per ROI

input: raw excel phenotype files after inital genetic / health etc QC
output: QC'ed excel pheno type files, sumstats, overview plots and final mask/eid list
'''
#import stuff
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os 
from statsmodels import robust
import csv

#------------------------------------------------------------------
# 1. PHENO QC
#------------------------------------------------------------------

# Set up paths
#change this depending on sample used
sample ='replication'
path_freesurfer='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/preprocessing' 
if not os.path.exists((path_freesurfer.format(sample)+ '/update')):
    os.mkdir((path_freesurfer.format(sample)+ '/update'))
#import the files
global_values = pd.read_excel ('{}/ukb43760_{}_globalValues.xlsx'.format(path_freesurfer.format(sample),sample), usecols = 'C:G')


#QC GLOBAL VALUES by median 5* MAD 
lower_bound =[]
upper_bound =[]
for column in global_values.iloc[:,0:]:
    lower = global_values[column].median() - 5*robust.mad(global_values[column])
    upper= global_values[column].median() + 5*robust.mad(global_values[column])
    lower_bound.append(lower)
    upper_bound.append(upper)
#empty dataframe
global_values_qc = pd.DataFrame()   
#check each column for outliers, make qc'ed dataframe
for ind, column in enumerate(global_values.iloc[:,0:].columns):
    global_values_qc[column] = global_values[column].mask((global_values[column] < lower_bound[ind]) | (global_values[column] > upper_bound[ind]))
# as global outliers get deleted, make a mask and save it 
global_nan = global_values_qc.notna()
mask_global=global_nan.all(axis=1)
mask_global.sum()
mask_global.to_csv('{}/update/mask_global_outlier_ukb43760.dat'.format(path_freesurfer.format(sample)), header=False, index =None)

# apply mask to your dataframe, we just want values there
global_values_qc = global_values_qc[mask_global].copy()
#get sum stats, add percentage ( note they are calculated from the TOTAL N before we deleted global outliers)
global_values_sumstats = global_values_qc.describe()
global_values_sumstats.loc['percentage',:] = global_values_sumstats.loc['count',:] / len(global_values.index)
#save stuff
global_values_sumstats.to_csv('{}/update/globalValues_QCed_sumstats.csv'.format(path_freesurfer.format(sample)))

# make violinplot with Seaborn - this here has to happen with subplots as we are plotting different metrics thus have different Y axis!
if not os.path.exists((path_freesurfer.format(sample) + '/update/plots')):
    os.mkdir((path_freesurfer.format(sample) + '/update/plots'))
f, axes = plt.subplots(1, 2, figsize=(10, 7), sharex=False)
sns.despine(left=True)
# Plot Thickness
sns.violinplot(y = 'value', x = 'variable' , data = pd.melt(global_values_qc.iloc[:,0:3:2]), palette = 'BuGn_r', ax=axes[0])
# Plot Surface area
sns.violinplot(y = 'value', x = 'variable' , data = pd.melt(global_values_qc.iloc[:,1:4:2]), palette = 'Blues', ax=axes[1])
plt.tight_layout()
axes[1].set_title("Thickness",fontsize=12)
axes[1].set_xlabel("ROI",fontsize=12)
axes[1].set_ylabel('thickness (mm)',fontsize=12)
axes[0].set_title("Surface area",fontsize=12)
axes[0].set_xlabel("ROI",fontsize=12)
axes[0].set_ylabel ('Area ($mm^2$)',fontsize=12)
f.tight_layout()
f.savefig('{}/update/plots/ukb43760_globalValues_QCed.png'.format(path_freesurfer.format(sample)),bbox_inches="tight")

#---------------------------------------------------------------------------------------------------------------------------
#now we start on the regional files derived via the DK atlas 
#import the files
thickness_left= pd.read_excel ('{}/ukb43760_{}_thickness.xlsx'.format(path_freesurfer.format(sample),sample), sheet_name = 'thickness_left', usecols = 'C:AI')
thickness_right= pd.read_excel ('{}/ukb43760_{}_thickness.xlsx'.format(path_freesurfer.format(sample),sample), sheet_name = 'thickness_right', usecols = 'C:AI')
surface_left = pd.read_excel ('{}/ukb43760_{}_surface.xlsx'.format(path_freesurfer.format(sample),sample), sheet_name = 'surface_left', usecols = 'C:AI')
surface_right = pd.read_excel ('{}/ukb43760_{}_surface.xlsx'.format(path_freesurfer.format(sample),sample), sheet_name = 'surface_right', usecols = 'C:AI')


lofdfs = [thickness_left, thickness_right, surface_left,surface_right]
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
        qc_global = qc[mask_global].copy()
        rep_pheno_qc[key]=qc_global.copy()

thickness_left_qc=rep_pheno_qc[0]
thickness_right_qc=rep_pheno_qc[1]
surface_left_qc=rep_pheno_qc[2]
surface_right_qc=rep_pheno_qc[3]

#-----------------------------------------------------
# quirk of pandas - can't work on columns in a global way otherwise so we save them to a list 
sLE_col =surface_left_qc.columns.values.tolist()
sRE_col =surface_right_qc.columns.values.tolist()
tLE_col =thickness_left_qc.columns.values.tolist()
tRE_col =thickness_right_qc.columns.values.tolist()

#and replace them with numbers
surface_left_qc.columns=np.arange(len(surface_left_qc.columns))
surface_right_qc.columns=np.arange(len(surface_right_qc.columns))
thickness_left_qc.columns=np.arange(len(thickness_left_qc.columns))
thickness_right_qc.columns=np.arange(len(thickness_right_qc.columns))

#-----------------------------------------------------------------
#now we replace outliers with NAN by checking each dataframe frame and get sumstats
thicknessLE_masked = thickness_left_qc.mask(surface_right_qc.isna() | surface_left_qc.isna() | thickness_right_qc.isna()).copy()
thicknessLEsumstats = thicknessLE_masked.describe() 
thicknessLEsumstats.loc['percentange',:] = thicknessLEsumstats.loc['count',:] / len(thickness_left.index)

thicknessRE_masked = thickness_right_qc.mask(surface_left_qc.isna() | surface_right_qc.isna() | thickness_left_qc.isna()).copy()
thicknessREsumstats = thicknessRE_masked.describe() 
thicknessREsumstats.loc['percentange',:] = thicknessREsumstats.loc['count',:] / len(thickness_right.index)

surfaceLE_masked = surface_left_qc.mask(surface_right_qc.isna() | thickness_left_qc.isna() | thickness_right_qc.isna()).copy()
surfaceLEsumstats = surfaceLE_masked.describe() 
surfaceLEsumstats.loc['percentange',:] = surfaceLEsumstats.loc['count',:] / len(surface_left.index)

surfaceRE_masked = surface_right_qc.mask(surface_left_qc.isna()| thickness_left_qc.isna() | thickness_right_qc.isna()).copy()
surfaceREsumstats = surfaceRE_masked.describe() 
surfaceREsumstats.loc['percentange',:] = surfaceREsumstats.loc['count',:] / len(surface_right.index)

#now name the columns again
surfaceLE_masked.columns =sLE_col
surfaceRE_masked.columns =sRE_col
thicknessRE_masked.columns =tRE_col
thicknessLE_masked.columns =tLE_col
thicknessLEsumstats.columns=tLE_col
thicknessREsumstats.columns=tRE_col
surfaceLEsumstats.columns=sLE_col
surfaceREsumstats.columns=sRE_col


#get EID list as it's needed later on for sorting 
eid_enigma = pd.read_csv('{}/eid_enigmaEvol_{}_ukb43760.dat'.format(path_freesurfer.format(sample),sample), header =None, names =['eid'],dtype=str)
eid_enigma = eid_enigma[mask_global].copy()
#save the final EID list for variant QC / subsetting etc
eid_enigma.to_csv('{}/update/eid_afterQC_{}_ukb43760.dat'.format(path_freesurfer.format(sample),sample),index=None,header=None)
#get the overall ID list to create final mask
eid = pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/current_release/ukb43760/ukb43760.csv', sep =",", usecols =[0])
mask_QCed = eid['eid'].isin(eid_enigma['eid'])
#let's print out the N just as a sanity check
mask_QCed.sum()
#save the final mask - this combines neuro healthy, imaging, pheno QC and genetic sample QC
mask_QCed.to_csv('{}/update/mask_afterQC_{}_ukb43760.dat'.format(path_freesurfer.format(sample),sample),index=None)
thickness_left_afterQC = pd.concat([eid_enigma,thicknessLE_masked],axis= 1)
thickness_right_afterQC = pd.concat([eid_enigma,thicknessRE_masked],axis= 1)
surface_left_afterQC = pd.concat([eid_enigma,surfaceLE_masked],axis= 1)
surface_right_afterQC = pd.concat([eid_enigma,surfaceRE_masked],axis= 1)
global_afterQC=pd.concat([eid_enigma,global_values_qc],axis=1)


#save the QC'ed files
writer = pd.ExcelWriter('{}/update/ukb43760_regionalDK_QCed_{}_sumstats.xlsx'.format(path_freesurfer.format(sample),sample), engine='xlsxwriter')
surfaceLEsumstats.to_excel(writer, sheet_name ='surface_LE')
surfaceREsumstats.to_excel(writer,sheet_name='surface_RE')
thicknessLEsumstats.to_excel(writer, sheet_name='thickness_LE')
thicknessREsumstats.to_excel(writer,sheet_name='thickness_RE')
writer.save()


writer = pd.ExcelWriter('{}/update/ukb43760_regionalDK_QCed_{}_values.xlsx'.format(path_freesurfer.format(sample),sample), engine='xlsxwriter')
surface_left_afterQC.to_excel(writer, sheet_name ='surface_LE')
surface_right_afterQC.to_excel(writer,sheet_name='surface_RE')
thickness_left_afterQC.to_excel(writer, sheet_name='thickness_LE')
thickness_right_afterQC.to_excel(writer,sheet_name='thickness_RE')
writer.save()

global_afterQC.to_excel('{}/update/ukb43760_globalValues_QCed_{}.xlsx'.format(path_freesurfer.format(sample),sample),sheet_name='global_values')

#---------------------------------------------------------------------------
# PLOT all regional metrics after QC as well 
#--------------------------------------------------------------------------
# make violinplot with Seaborn
# loop over it to make it more readable 

plot_path='{}/update/plots/ukb43760_{}_QCed.png'
metrics={'thicknessLE':thicknessLE_masked,'thicknessRE': thicknessRE_masked}

for key, value in metrics.items():
    plt.figure(figsize=(20,5))
    bplot=sns.violinplot(y='value', x='variable', 
                     data=pd.melt(value), 
                     width=0.5,
                     palette="GnBu_d")
    bplot.set_xticklabels(bplot.get_xticklabels(), rotation=45, horizontalalignment='right')
    sns.despine()
    bplot.set_xlabel("ROI",fontsize=12)
    bplot.set_ylabel("Thickness (mm)",fontsize=12)
    bplot.tick_params(labelsize=12)
    bplot.autoscale()
    fig = bplot.get_figure()
    fig.savefig(plot_path.format(path_freesurfer.format(sample),key),bbox_inches="tight")

    
metrics={'surfaceLE': surfaceLE_masked,'surfaceRE': surfaceRE_masked}
for key, value in metrics.items():
# make violinplot with Seaborn
    plt.figure(figsize=(20,5))
    bplot=sns.violinplot(y='value', x='variable', 
                     data=pd.melt(value), 
                     width=0.5,
                     palette="BuGn_r")
    bplot.set_xticklabels(bplot.get_xticklabels(), rotation=45, horizontalalignment='right')
    sns.despine()
    bplot.set_xlabel("ROI",fontsize=12)
    bplot.set_ylabel("Surface area ($mm^2$)",fontsize=12)
    bplot.tick_params(labelsize=12)
    bplot.autoscale()
    fig = bplot.get_figure()
    fig.savefig(plot_path.format(path_freesurfer.format(sample),key),bbox_inches="tight")

#---------------------------------------------------------------------------
# 2 Prepare Covariates 
#--------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
# get the global values  --> is here now global_values_qc
# get the other covariates used in the analysis 
confound_values =pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/current_release/ukb43760/ukb43760.csv', sep =",", usecols =['eid','31-0.0','21003-2.0','22009-0.1','22009-0.2','22009-0.3','22009-0.4','22009-0.5','22009-0.6','22009-0.7','22009-0.8', '22009-0.9', '22009-0.10', '25756-2.0', '25757-2.0', '25758-2.0', '25734-2.0','25735-2.0','54-2.0','22000-0.0'])
# sort the columns again as pandas rearranges them... 
confound_values = confound_values[['eid','31-0.0','21003-2.0','22009-0.1','22009-0.2','22009-0.3','22009-0.4','22009-0.5','22009-0.6','22009-0.7','22009-0.8', '22009-0.9', '22009-0.10', '25756-2.0', '25757-2.0', '25758-2.0', '25734-2.0','25735-2.0','54-2.0','22000-0.0']]
confound_values.columns = ['eid','sex','age_at_scan','PC1','PC2','PC3','PC4','PC5', 'PC6','PC7', 'PC8','PC09','PC10','xposition','yposition','zposition','t1SNR','t1CNR','assessmentcentre','genotypearray']
#change sex value from 0 female 1 male to female = 2 male =1 
confound_values.loc[confound_values['sex'] == 0, 'sex'] = 2
#make genotype batch into a binary array variables
confound_values.loc[confound_values['genotypearray'] < 0, 'genotypearray'] = 0
confound_values.loc[confound_values['genotypearray'] > 0, 'genotypearray'] = 1
# get mean age in sample
age_mean =confound_values['age_at_scan'].mean()
# get mean centred age
confound_values.loc[:,'ageCSq'] = confound_values.loc[:,'age_at_scan'].apply(lambda x: np.square(x - age_mean))
# get sex age interaction
confound_values.loc[:,'age_sex'] = confound_values.loc[:,'age_at_scan'] * confound_values.loc[:,'sex']
#get sex mean centred age interaction
confound_values.loc[:,'ageCSq_sex'] = confound_values.loc[:,'ageCSq'] * confound_values.loc[:,'sex']
#let's factorize the categorical variables 
cat_variables=['sex','assessmentcentre','genotypearray']
for col in cat_variables:
    confound_values[col] = pd.factorize(confound_values[col])[0]
#as we are manipulating the dataframe, make a copy
#also as we dummy code to this after subsetting the inital array as we otherwise end up with more dummy variables than needed
confound_sample = confound_values[mask_QCed].copy()
confound_sample=pd.get_dummies(confound_sample, columns=['assessmentcentre'],drop_first =True)
confound_sample['eid']=confound_sample.eid.astype(str)
#---------------------------------------------------------------------------------------------------------------
#save it
#---------------------------------------------------------------------------------------------------------------
#WITHOUT GLOBAL FS metrics
confound_sample.to_excel('{}/update/ukb43760_covariats_withoutGlob_{}.xlsx'.format(path_freesurfer.format(sample),sample),index = None)
if sample == 'replication':
    # REPLICATIION create global covariates
    global_covariates_rep=pd.DataFrame()
    global_values_qc.reset_index(inplace=True,drop=True)
    global_covariates_rep['global_surface'] = global_values_qc['26721-2.0'] + global_values_qc['26822-2.0']
    global_covariates_rep['global_thickness'] =0.5*( global_values_qc['26755-2.0'] + global_values_qc['26856-2.0'])
    global_covariates_rep['global_thickness']=global_covariates_rep['global_thickness'].round(4)
    confound_sample.reset_index(inplace =True,drop=True)
    confound_sample_thickness = pd.concat([confound_sample, global_covariates_rep['global_thickness']],axis=1)
    confound_sample_surface = pd.concat([confound_sample, global_covariates_rep['global_surface']],axis=1)
    confound_sample_thickness.to_excel('{}/update/ukb43760_covariats_thickness_{}.xlsx'.format(path_freesurfer.format(sample),sample),index = None)
    confound_sample_surface.to_excel('{}/update/ukb43760_covariats_surface_{}.xlsx'.format(path_freesurfer.format(sample),sample),index = None)
    #---------------------------------------------------------------------------------------------------------------#save this
#EU - hemisphere specific
#---------------------------------------------------------------------------------------------------------------#save this
if sample == 'european':
    confound_sample.reset_index(inplace =True,drop=True)
    global_values_qc.reset_index(inplace=True,drop=True)
    confound_sample_thickness_le = pd.concat([confound_sample, global_values_qc['26755-2.0']],axis=1)
    confound_sample_thickness_re = pd.concat([confound_sample, global_values_qc['26856-2.0']],axis=1)
    confound_sample_surface_le = pd.concat([confound_sample, global_values_qc['26721-2.0']],axis=1)
    confound_sample_surface_re = pd.concat([confound_sample, global_values_qc['26822-2.0']],axis=1)
    #save this
    confound_sample_thickness_le.to_excel('{}/update/ukb43760\-covariats_thickness_le_{}.xlsx'.format(path_freesurfer.format(sample),sample),index = None)
    confound_sample_thickness_re.to_excel('{}/update/ukb43760_covariats_thickness_re_{}.xlsx'.format(path_freesurfer.format(sample),sample),index = None)
    confound_sample_surface_le.to_excel('{}/update/ukb43760_covariats_surface_le_{}.xlsx'.format(path_freesurfer.format(sample),sample),index = None)
    confound_sample_surface_re.to_excel('{}/update/ukb43760_covariats_surface_re_{}.xlsx'.format(path_freesurfer.format(sample),sample),index = None)


#----------------------------------------------------------------------------------
# 3 BEGENIE prep
#------------------------------------------------------------------------------------

# PHENOTYPES
path_bgenie='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/bgenieInput' 

if not os.path.exists(path_bgenie.format(sample)):
    os.mkdir(path_bgenie.format(sample))
#read the sample file, strip first row ( 0 0 0) and only import one column
bgen_list = pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/snpQC/subsetting/ukb43760_enigmaEvo_{}_chr2.sample'.format(sample,sample),delim_whitespace=True, usecols=[0],skiprows=[1],dtype=str)
# read the roi list - needed later to give better column names (currently UKB IDPs)
roi_name = pd.read_excel('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/QC_5MAD/ROI_table.xlsx', header =None)
#--------------------------------
#--------------------------------
global_afterQC.set_index('eid', inplace=True)
# sort by the sample file
global_pheno_sorted = global_afterQC.loc[bgen_list['ID_1']]

regional_fs = {'surface_le':surface_left_afterQC, 'surface_re':surface_right_afterQC, 'thickness_le':thickness_left_afterQC,'thickness_re':thickness_right_afterQC}
for key,value in regional_fs.items():
   regional_fs[key].set_index('eid', inplace=True)
   regional_fs[key]=regional_fs[key].loc[bgen_list['ID_1']]    
   regional_fs[key].columns=np.arange(len(regional_fs[key].columns))    

#----------------------------
# REPLICATION sample
#--------------------------------
if sample == 'replication':
#for replication we neet total SA and mean TH
    global_replication=pd.DataFrame()
    global_replication['global_surface'] = (global_pheno_sorted['26721-2.0'] + global_pheno_sorted['26822-2.0'])
    global_replication['global_thickness'] =0.5*( global_pheno_sorted['26755-2.0'] + global_pheno_sorted['26856-2.0'])
    global_replication= global_replication.round(decimals=4)
    #save as a white space delimitered table 
    global_replication.to_csv('{}/ukb43760_global_{}.table'.format(path_bgenie.format(sample),sample), index = None, sep=' ', quoting = csv.QUOTE_NONE)
    
#also average and set up the regional measures     
    regional_replication = dict.fromkeys(['thickness','surface'])
    regional_replication['surface']= 0.5 *(regional_fs['surface_le'] +regional_fs['surface_re'] )
    regional_replication['thickness'] = 0.5 *(regional_fs['thickness_le'] +regional_fs['thickness_le'])
    for key, value in regional_replication.items():
        regional_replication[key].fillna(-999, inplace=True)
        regional_replication[key].columns = roi_name.iloc[:,0]
        regional_replication[key]=value.round(decimals=4)
        #save as a white space delimitered table 
        regional_replication[key].to_csv('{}/ukb43760_regional_{}_{}.table'.format(path_bgenie.format(sample),key,sample), index = None, sep=' ', quoting = csv.QUOTE_NONE)
        
#--------------------------
# LEFT / RIGHT EU sample
#--------------------------
if sample == 'european':    
    global_european= global_pheno_sorted .round(decimals=4)
    #save as a white space delimitered table 
    global_european.to_csv('{}/ukb43760_global_hemi_{}.table'.format(path_bgenie.format(sample),sample), index = None, sep=' ', quoting = csv.QUOTE_NONE) 
    
    for key,value in regional_fs.items():
        regional_fs[key].fillna(-999, inplace=True)
        regional_fs[key]=regional_fs[key].round(decimals=4)
        roi_name_new=roi_name[0] +'_'+key[-2:]
        regional_fs[key].columns = roi_name_new
        regional_fs[key].to_csv('{}/ukb43760_regional_{}_{}.table'.format(path_bgenie.format(sample),key,sample), index = None, sep=' ', quoting = csv.QUOTE_NONE)


#--------------------------------
# COVARIATES
#-------------------------------
output_name='{}/ukb43760_covariates_{}_{}.table'

if sample == 'replication':
    confound_replication =[confound_sample, confound_sample_surface,confound_sample_thickness] 
    for i,metric in enumerate(['noGlob','surface', 'thickness']):
        confound_replication[i].set_index('eid', inplace=True)
        confound_replication[i]=confound_replication[i].loc[bgen_list['ID_1']]
        confound_replication[i]= confound_replication[i].round(decimals=4)
        confound_replication[i].to_csv(output_name.format(path_bgenie.format(sample),metric,sample), index = None, sep=' ', quoting = csv.QUOTE_NONE)

if sample == 'european':
    confound_hemi =[confound_sample, confound_sample_surface_le,confound_sample_surface_re,confound_sample_thickness_le,confound_sample_thickness_re]        
    for i,metric in enumerate(['noGlob','surface_le', 'surface_re','thickness_le', 'thickness_re']):
        confound_hemi[i].set_index('eid', inplace=True)
        confound_hemi[i] =confound_hemi[i].loc[bgen_list['ID_1']]
        confound_hemi[i]= confound_hemi[i].round(decimals=4)
        confound_hemi[i].to_csv(output_name.format(path_bgenie.format(sample),metric,sample), index = None, sep=' ', quoting = csv.QUOTE_NONE)       
            
  
#--------------------------
# SNPs to keep - original files created by variantQC procedure
#--------------------------
file_name = '/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/snpQC/variantQC/ukb43760_enigmaEvo_{}_chr{}.snpstats_mfi_hrc.snps2keep'
output_name = '{}/snps2keep/ukb43760_enigmaEvo_{}_snps2keep_chr{}.table'

if not os.path.exists((path_bgenie.format(sample)+ '/snps2keep')):
    os.mkdir((path_bgenie.format(sample)+ '/snps2keep'))

for i in range(1, 23):
    current_data=pd.read_csv(file_name.format(sample,sample,i),sep='\s+',usecols =['RS_ID.UKB'])
    current_data.to_csv(output_name.format(path_bgenie.format(sample),sample,i), index = None, sep=' ', quoting = csv.QUOTE_NONE)