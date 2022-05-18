#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 11:50:46 2021
Modified May 2022
@author: B. Molz
"""
# This script uses either the LDSC genetic correlation summary .txt file to calculate corrleation matrices
# raw matrices + heatmaps are save

#import stuff
import pandas as pd
import csv
import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np
#-------------------------------------------------------------------------------------------------------------------
#first work on genetic correlations
# input file name
file_name ='/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/le_re/ldsc/gencor_phenoSPD/gencor_surface_phenoSPD_allin.txt'
# also need filenames to merge the chunks later without issues
trait_name =pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/le_re/ldsc/gencor_phenoSPD/trait_names.txt',sep='\s+',usecols=['trait'])

#depending on your current trait names, you might have to modify the list slighty to fit the order you imported the traits
# create new metric
trait_name['hemi'] = np.where(trait_name['trait'].str.contains('le',case=False), 'le', 're')
# now we delete the inital hemisphere info in the trait as our current traits have this info after the ROI info
trait_name['trait']=trait_name['trait'].str.replace('_le','').str.replace('_re','')
# create the new name
trait_name['new'] =trait_name['trait'] + '_' + trait_name['hemi']
#sort this alphabetically - this should now be identical with the imported order of traits
trait_name=trait_name.sort_values('new',ascending=True)

#empty list to hold our chuncks
df_list =[]

#loop over the gen cor file, include chuncks into the list
for chunk in pd.read_csv(file_name, sep='\s+', chunksize=68, usecols=['rg']):
    df_list.append(chunk)
    
#posibbly index to different trait_name.iloc  
#to concat chucks, reset the index, and give them individual column names
for i in range(0,68):
    df_list[i].reset_index(drop=True,inplace=True)
    df_list[i].rename(columns={'rg': trait_name.iloc[i,0]},inplace=True)
    
#concat 
matrix_gen = pd.concat(df_list,axis=1)
#also include rownames
matrix_gen.index=trait_name
#plot and save the heatmap
f, ax=plt.subplots(figsize=(20,20))
sn.heatmap(matrix_gen, fmt="g", center=0, cmap=sn.diverging_palette(245, 5, as_cmap=True),vmin= -1, vmax=1,cbar_kws={ "ticks":[-1,-0.5,0,0.5,1],"shrink": .30})
bottom, top = ax.get_ylim()
ax.set_ylim(bottom + 0.5, top - 0.5)
f.savefig('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/le_re/ldsc/genCor_DTI_heatmap.png',bbox_inches="tight")
#save the matrix without any headers/index labels 
matrix_gen.to_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/le_re/ldsc/genCor_surfaceHemi.txt', sep=' ', quoting = csv.QUOTE_NONE,header=False, index=False)
matrix_gen.to_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/european/le_re/ldsc/gencor_surfaceHemi_withLabels.txt', sep=' ', quoting = csv.QUOTE_NONE)

