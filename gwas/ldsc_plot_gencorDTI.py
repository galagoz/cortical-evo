'''
This script uses the genetic correlation table to plot r(g) values against the specific traits, here surface area and thickness 
Input: table with r(g), se, p values etc
Output: graph 
'''
#import stuff
import pandas as pd 
import matplotlib.pyplot as plt
#---------------------------------------------------------------------------------------------------------
#set up the working sample
sample='european'
#read the summary gen cor file
gencor = pd.read_csv('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/DTI/LDSC/gencor/gencor_ukbBig40_DTI.txt'.format(sample), delim_whitespace=True)
# get rid of gencors without global FS correction 
#replace stuff from trait we don'want in the graph
gencor['trait']=gencor['trait'].str.replace('gencor_ukb43760_DTI_','')
#format the pValue
gencor.p = gencor.p.map(lambda x: '{:.2E}'.format(x))
gencor['traitP'] = gencor['trait'] + '  p=['+ gencor['p'].astype(str) + ' ]'


#plot
fig,ax1 = plt.subplots(figsize=(5,20))
ax1.errorbar(gencor.rg, gencor.traitP, xerr=gencor.se, fmt='o', color='green',
             ecolor='darkgrey', elinewidth=3, capsize=0, label = 'DTI');
ax1.axvline(x=1,ymin=0,ymax=1, color='red', ls=':', lw=2)         
ax1.grid('on')
ax1.set_xlabel('rg', fontsize =12)
ax1.tick_params(axis='both', which='minor', labelsize=12)
ax1.margins(y=0.03)
plt.xlim(0.5, 1.2)

fig.savefig('/data/workspaces/lag/workspaces/lg-ukbiobank/projects/enigma_evol/final/{}/DTI/LDSC/gencor.png'.format(sample),bbox_inches="tight")