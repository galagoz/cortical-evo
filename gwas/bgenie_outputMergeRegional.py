#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 13:34:52 2021

@author: Barbara Molz
"""
'''
This script loads the output created by BGNENIE, depending on phenotypes (and possible differences in covariates), merges all chromosomes, converts the -log10p value back to normal p values and saves it as a CSV
Input: path to input files, path where output should be saved, phenotypes
Output: CSV folder per pheno with all chromosomes
'''
import pandas as pd
import csv

sample='european' 
# create the input and output path variables
#workingpath
workingpath='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/bgenieOutput/{}/{}'
#file_name ='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/averaged/bgenie/bgenie_ukb43760_globalValues_chr{}_{}.out.gz'
filename_rep='{}/bgenie_ukb43760_regional_{}_{}_chr{}_{}.out.gz'
filename_eu='{}/bgenie_ukb43760_regional_{}_{}_{}_chr{}_{}.out.gz'

#output_name='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/averaged/sumstats/sumstats_ukb43760_{}_{}_allChr.gz'
outputName='{}/sumstats_ukb43760_regional_{}_{}_{}_{}_allChr.gz' 
n_count= pd.read_excel('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/preprocessing/update/ukb43760_regionalDK_QCed_{}_sumstats.xlsx'.format(sample,sample),nrows =1, usecols ='B:AH')
roi_name = pd.read_excel('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/enigma_evo/QC_5MAD/ROI_table.xlsx', header =None)
n_count.columns=roi_name[0]

#----------------------------
#empty list to store the seperate dataframes (per chromosom)
# get the N for each phenotype
df_list =[]
df_list2=[]

if sample == 'european':
    for pheno in (['surface','thickness']):
        for hemi in (['le', 're']):
            for covar in (['withGlob', 'noGlob']):
                for i in range(0, 22):
                    #phenotype, just load what's true for all sub ROIs: Chr, position, rsID, allele 1+2, info, af
                    df_list.append(pd.read_csv(filename_eu.format(workingpath.format(sample,pheno,covar),pheno,hemi,covar,i+1,sample),sep='\s+', usecols =[*range(0,7)]))
                df_list = pd.concat(df_list)
                #create the marker
                df_list['marker']=df_list['chr'].astype(str) +':'+ df_list['pos'].astype(str)
                #rename columns
                df_list.columns = ['CHR','SNP', 'BP', 'A2', 'A1', 'FREQ1','INFO','MARKER']
                for roi in range(7,136,4):
                    for chrom in range(0, 22):
                        #now get the phenotypes, one at at a time, we want: beta, SE and -log10p
                        df_list2.append(pd.read_csv(filename_eu.format(workingpath.format(sample,pheno,covar),pheno,hemi,covar,chrom+1,sample),sep='\s+', usecols =[roi, roi+1,roi+3]))
                    df_list2 = pd.concat(df_list2)
                    #convert the column names for the -log10 pvalues
                    col_names =df_list2.columns[2]
                    #convert -log10p values to normal p values again - easier to handle later on
                    df_list2[col_names]=10**(-df_list2[col_names])
                    #not neccessary here but convenient to use for output file name
                    col_index=col_names.replace('-log10p','')
                    #we need the correct N - this differs per sample 
                    df_list2.insert(3,"N",n_count.at[0,col_index[:-3]], allow_duplicates=False)
                    #rename the columns
                    df_list2.columns = ['BETA', 'SE', 'P', 'N']
                    #concatenate this with the main dataframe 
                    df_final= pd.concat([df_list,df_list2],axis=1)
                    #reorder columns
                    df_final = df_final[['SNP', 'A1', 'A2', 'FREQ1', 'BETA', 'SE', 'P', 'N', 'MARKER', 'CHR',  'BP','INFO']]
                    #save the dataframe and empty the list again as a zipped file using the DK name
                    df_final.to_csv(outputName.format(workingpath.format(sample,pheno,covar),pheno,col_index,covar,sample),index = None,sep=' ', quoting = csv.QUOTE_NONE, compression='gzip')
                    #empty the phenotpe dataframe and start again
                    df_list2=[]
                    del df_final
                df_list=[]
else:
    for pheno in (['surface','thickness']):
        for covar in (['withGlob', 'noGlob']):
            for i in range(0, 22):
                #phenotype, just load what's true for all sub ROIs: Chr, position, rsID, allele 1+2, info, af
                df_list.append(pd.read_csv(filename_rep.format(workingpath.format(sample,pheno,covar),pheno,covar,i+1,sample),sep='\s+', usecols =[*range(0,7)]))
            df_list = pd.concat(df_list)
            #create the marker
            df_list['marker']=df_list['chr'].astype(str) +':'+ df_list['pos'].astype(str)
            #rename columns
            df_list.columns = ['CHR','SNP', 'BP', 'A2', 'A1', 'FREQ1','INFO','MARKER']
            for roi in range(7,136,4):
                for chrom in range(0, 22):
                    #now get the phenotypes, one at at a time, we want: beta, SE and -log10p
                    df_list2.append(pd.read_csv(filename_rep.format(workingpath.format(sample,pheno,covar),pheno,covar,chrom+1,sample),sep='\s+', usecols =[roi, roi+1,roi+3]))
                df_list2 = pd.concat(df_list2)
                #convertt the column names for the -log10 pvalues
                col_names =df_list2.columns[2]
                #convert -log10p values to normal p values again - easier to handle later on
                df_list2[col_names]=10**(-df_list2[col_names])
                #not neccessary here but convenient to use for output file name
                col_index=col_names.replace('-log10p','')
                #we need the correct N - this differs per sample 
                df_list2.insert(3,"N",n_count.at[0,col_index], allow_duplicates=False)
                #rename the columns
                df_list2.columns = ['BETA', 'SE', 'P', 'N']
                #concatante this with the main dataframe 
                df_final= pd.concat([df_list,df_list2],axis=1)
                #reorder columns
                df_final = df_final[['SNP', 'A1', 'A2', 'FREQ1', 'BETA', 'SE', 'P', 'N', 'MARKER', 'CHR',  'BP','INFO']]
                #save the dataframe and empty the list again as a zipped file using the DK name
                df_final.to_csv(outputName.format(workingpath.format(sample,pheno,covar),pheno,col_index,covar,sample),index = None,sep=' ', quoting = csv.QUOTE_NONE, compression='gzip')
                #empty the phenotpe dataframe and start again
                df_list2=[]
                del df_final
            df_list=[]