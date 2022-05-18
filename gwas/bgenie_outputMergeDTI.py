#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created Nov 2021
@author: Barbara Molz
"""
'''
This script loads the output created by BGNENIE, depending on phenotypes (and possible differences in covariates), merges all chromosomes, converts the -log10p value back to normal p values and saves it as a CSV
Input: path to input files, path where output should be saved, phenotypes
Output: CSV folder per pheno with all chromosomes
'''
import pandas as pd
import csv
import os

sample='europeanDTI' 
# create the input and output path variables
pathDTI_gwas='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/DTI/GWAS' 
#file_name ='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/averaged/bgenie/bgenie_ukb43760_globalValues_chr{}_{}.out.gz'
file_name='{}/bgenie_ukb43760_{}_chr{}.out.gz'
#output_name='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/averaged/sumstats/sumstats_ukb43760_{}_{}_allChr.gz'
pathDTI_sumstats='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/DTI/GWAS/sumstats' 
if not os.path.exists(pathDTI_sumstats):
    os.mkdir(pathDTI_sumstats)
output_name='{}/sumstats_ukb43760_{}_{}_allChr.gz' 
n_count=pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/DTI/ukb43760_dtiTBSS_QCed_sumstats.csv'.format(sample),nrows=1)
UKB_idps = pd.read_html('/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/current_release/ukb43760/ukb43760.html',header=1)
UKB_idps = pd.read_html('/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/current_release/ukb43760/ukb43760.html')
UKB_idps = UKB_idps[1]
#----------------------------
#empty list to store the seperate dataframes (per chromosom)
# get the N for each phenotype
df_list =[]
df_list2=[]
# lets only import what we need as files are already massive
for i in range(0, 22):
    #per overall phenotype, just load what's true for all sub ROIs: Chr, position, rsID, allele 1+2, info, af
    df_list.append(pd.read_csv(file_name.format(pathDTI_gwas,sample,i+1,sample),sep='\s+', usecols =[*range(0,7)]))
df_list = pd.concat(df_list)
#create the marker
df_list['marker']=df_list['chr'].astype(str) +':'+ df_list['pos'].astype(str)
#rename columns
df_list.columns = ['CHR','SNP', 'BP', 'A2', 'A1', 'FREQ1','INFO','MARKER']
# as we don't want to load the contant columns per pheno again, make another loop
for roi in range(7,196,4):
    for l in range(0, 22):
        #now get the phenotypes, one at at a time, we want: beta, SE and -log10p
        df_list2.append(pd.read_csv(file_name.format(pathDTI_gwas,sample,l+1),sep='\s+', usecols =[roi, roi+1,roi+3]))
        #concatenate the list into one single data frame
    df_list2 = pd.concat(df_list2)
    #get the column names for the -log10 pvalues
    col_names =df_list2.columns[2]
    #convert -log10p values to normal p values again - easier to handle later on
    df_list2[col_names]=10**(-df_list2[col_names])
    #not neccessary here but convenient to use for output file name
    col_index=col_names.replace('-log10p','')
    traitName=UKB_idps[UKB_idps['UDI']==col_index]['Description'].values[0]
    traitName=traitName.replace("Mean FA in ","").replace(" on FA skeleton","").replace(" ","_").replace('(','').replace(')','')
    #we need the correct N - this differs per sample 
    df_list2.insert(3,"N",int(n_count.at[0,col_index]), allow_duplicates=False)
    #rename the columns
    df_list2.columns = ['BETA', 'SE', 'P', 'N']
    #concatante this with the main dataframe 
    df_final= pd.concat([df_list,df_list2],axis=1)
    #reorder columns
    df_final = df_final[['SNP', 'A1', 'A2', 'FREQ1', 'BETA', 'SE', 'P', 'N', 'MARKER', 'CHR',  'BP','INFO']]
    #save the dataframe and empty the list again as a zipped file using the DK name
    df_final.to_csv(output_name.format(pathDTI_sumstats,sample,traitName),index = None,sep=' ', quoting = csv.QUOTE_NONE, compression='gzip')
    #empty the phenotpe dataframe and start again
    df_list2=[]
    del df_final
df_list=[]     