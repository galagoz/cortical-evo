'''
This script uses the genetic correlation table to plot r(g) values against the specific traits, here surface area and thickness 
Input: table with r(g), se, p values etc
Output: graph 
Created by B. Molz, May 2022
'''
#import stuff
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
#---------------------------------------------------------------------------------------------------------
#set up the working sample
sample='replication'
#read the summary gen cor file
gencor = pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/LDSCUpdate/gencor_E3vsRep/summary_gencor_E3.txt'.format(sample), delim_whitespace=True)
# get rid of gencors without global FS correction 
gencor_withGlob=gencor[~gencor.trait.str.contains('withoutGlob')].copy()
#add a column with metric - not realy needed but we delete metric later from trait name so might be handy if things go wrong
gencor_withGlob['metric'] = np.where(gencor_withGlob['trait'].str.contains('surf',case=False), 'surface', 'thickness')
#replace stuff from trait we don'want in the graph
gencor_withGlob['trait']=gencor_withGlob['trait'].str.replace('surface_','').str.replace('gencorE3_Rep_','').str.replace('_withGlob','')
#format the pValue
gencor_withGlob.p = gencor_withGlob.p.map(lambda x: '{:.2E}'.format(x))
gencor_withGlob['traitP'] = gencor_withGlob['trait'] + '  p=['+ gencor_withGlob['p'].astype(str) + ' ]'

#sort stuff
surface =gencor_withGlob[gencor_withGlob.metric.str.contains('surface')].copy()
thickness =gencor_withGlob[gencor_withGlob.metric.str.contains('thickness')].copy()

#plot
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(11,9))

fig.subplots_adjust(wspace=2.12)

fig,ax1 = plt.subplots(figsize=(5,20))
ax1.errorbar(surface.iloc[34:0:-1,1], surface.iloc[34:0:-1,12], xerr=surface.iloc[34:0:-1,2], fmt='o', color='green',
             ecolor='darkgrey', elinewidth=3, capsize=0, label = 'surface');
ax1.errorbar(surface.iloc[0,1], [surface.iloc[0,12]], xerr=surface.iloc[0,2], fmt='o', color='black',
             ecolor='green', elinewidth=4, capsize=0);
#ax2.errorbar(thickness.iloc[34:0:-1,1], thickness.iloc[34:0:-1,12], xerr=thickness.iloc[34:0:-1,2], fmt='o', color='red',
#             ecolor='lightgray', elinewidth=3, capsize=0, label = 'thickness');
#ax2.errorbar(thickness.iloc[0,1], [thickness.iloc[0,12]], xerr=thickness.iloc[0,2], fmt='o', color='black',
##             ecolor='red', elinewidth=4, capsize=0);
#ax1.axvline(x=1,ymin=0,ymax=1, color='red', ls=':', lw=2)         
#ax2.axvline(x=1,ymin=0,ymax=1, color='red', ls=':', lw=2)             
ax1.grid('on')
ax1.set_xlabel('rg', fontsize =12)
ax1.tick_params(axis='both', which='minor', labelsize=12)
ax1.margins(y=0.03)
#ax2.grid('on')
#ax2.set_xlabel('rg', fontsize =12)
#ax2.tick_params(axis='both', which='minor', labelsize=12)
#ax2.margins(y=0.02)
fig.savefig('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/replication/LDSCUpdate/gencor_E3vsRep/genCor_E3.png',bbox_inches="tight")