#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created Dec 2021
@author: Barbara Molz
"""
'''
This script loads the output created by BGNENIE, depending on phenotypes (and possible differences in covariates), merges all chromosomes, converts the -log10p value back to normal p values and saves it as a CSV
Input: path to input files, path where output should be saved, phenotypes
Output: CSV folder per pheno with all chromosomes
'''
import pandas as pd
import csv

sample='replication' 
# create the input and output path variables
#workingpath
workingpath='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/bgenieOutput/global/'
#file_name ='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/averaged/bgenie/bgenie_ukb43760_globalValues_chr{}_{}.out.gz'
filename='{}/bgenie_ukb43760_global_chr{}_{}.out.gz'
#output_name='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/averaged/sumstats/sumstats_ukb43760_{}_{}_allChr.gz'
output_name='{}/sumstats_ukb43760_global_{}_{}_allChr.gz' 
n_count=pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/preprocessing/update/ukb43760_globalValues_QCed_{}_sumstats.csv'.format(sample,sample),nrows=1, index_col=[0])
globalID=pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/global_IDs.csv',sep=';')
#----------------------------
#empty list to store the seperate dataframes (per chromosom)
# get the N for each phenotype
df_list =[]
df_list2=[]
# lets only import what we need as files are already massive
for i in range(0, 22):
    #per overall phenotype, just load what's true for all sub ROIs: Chr, position, rsID, allele 1+2, info, af
    df_list.append(pd.read_csv(filename.format(workingpath.format(sample),i+1,sample),sep='\s+', usecols =[*range(0,7)]))
df_list = pd.concat(df_list)
#create the marker
df_list['marker']=df_list['chr'].astype(str) +':'+ df_list['pos'].astype(str)
#rename columns
df_list.columns = ['CHR','SNP', 'BP', 'A2', 'A1', 'FREQ1','INFO','MARKER']
# as we don't want to load the contant columns per pheno again, make another loop

if sample == 'european':
    for roi in range(7,20,4):
        for chrom in range(0, 22):
            #now get the phenotypes, one at at a time, we want: beta, SE and -log10p
            df_list2.append(pd.read_csv(filename.format(workingpath.format(sample),chrom+1,sample),sep='\s+', usecols =[roi, roi+1,roi+3]))
        df_list2 = pd.concat(df_list2)
        #convertt the column names for the -log10 pvalues
        col_names =df_list2.columns[2]
        #convert -log10p values to normal p values again - easier to handle later on
        df_list2[col_names]=10**(-df_list2[col_names])
        #not neccessary here but convenient to use for output file name
        col_index=col_names.replace('-log10p','').replace('global_','')
        col_index=globalID[globalID['ID']==col_index]['name'].values[0]
        #we need the correct N - this differs per sample 
        df_list2.insert(3,"N",n_count.values[0][0], allow_duplicates=False)
        #rename the columns
        df_list2.columns = ['BETA', 'SE', 'P', 'N']
        #concatante this with the main dataframe 
        df_final= pd.concat([df_list,df_list2],axis=1)
        #reorder columns
        df_final = df_final[['SNP', 'A1', 'A2', 'FREQ1', 'BETA', 'SE', 'P', 'N', 'MARKER', 'CHR',  'BP','INFO']]
        #save the dataframe and empty the list again as a zipped file using the DK name
        df_final.to_csv(output_name.format(workingpath.format(sample),col_index,sample),index = None,sep=' ', quoting = csv.QUOTE_NONE, compression='gzip')
        #empty the phenotpe dataframe and start again
        df_list2=[]
        del df_final
else:
    for roi in range(7,12,4):
        for chrom in range(0, 22):
            #concatenate the list into one single data frame
            df_list2.append(pd.read_csv(filename.format(workingpath.format(sample),chrom+1,sample),sep='\s+', usecols =[roi, roi+1,roi+3]))
        df_list2 = pd.concat(df_list2)
        #convertt the column names for the -log10 pvalues
        col_names =df_list2.columns[2]
        #convert -log10p values to normal p values again - easier to handle later on
        df_list2[col_names]=10**(-df_list2[col_names])
        #not neccessary here but convenient to use for output file name
        col_index=col_names.replace('-log10p','').replace('global_','')
        #we need the correct N - this differs per sample 
        df_list2.insert(3,"N",n_count.values[0][0], allow_duplicates=False)
        #rename the columns
        df_list2.columns = ['BETA', 'SE', 'P', 'N']
        #concatante this with the main dataframe 
        df_final= pd.concat([df_list,df_list2],axis=1)
        #reorder columns
        df_final = df_final[['SNP', 'A1', 'A2', 'FREQ1', 'BETA', 'SE', 'P', 'N', 'MARKER', 'CHR',  'BP','INFO']]
        #save the dataframe and empty the list again as a zipped file using the DK name
        df_final.to_csv(output_name.format(workingpath.format(sample),col_index,sample),index = None,sep=' ', quoting = csv.QUOTE_NONE, compression='gzip')
        #empty the phenotpe dataframe and start again
        df_list2=[]
        del df_final
df_list=[]     