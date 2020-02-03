#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 10:53:23 2020

@author: meike
"""

'''
Density histograms of genome sizes. Makes one figure with 4 subplots: One all three merged and one per genus
'''

import os
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

def get_info(file, field):
    '''
    Give input file and column of interest. Returns list with all items of that field.
    '''
    item_l = []
    with open(file) as f:
        headers = f.readline().strip().split('\t')
        headers_inds = {name: i for i, name in enumerate(headers)}
            
        for line in f:
            a = line.strip().split('\t')
            for i, item in enumerate(a):
                if headers_inds[field] == i and item != '':
                    if item.isdigit()==1:
                        item_l.append(int(item))
                    else:
                        item_l.append(item)
    return item_l


    
    
path = os.getcwd()
p = Path(path)
#Get the gs sizes in Mb from all species
lactofile = os.path.join(p.parents[0], 'files', '06012020_lactococcus_database.tsv')
lac_gs = get_info(lactofile, 'genome_length')
lac_gs = [round(elem/1000000, 2) for elem in lac_gs]

streptofile = os.path.join(p.parents[0], 'files', '20012020_streptococcus_database.tsv')
strepto_gs = get_info(streptofile, 'genome_length')
strepto_gs = [round(elem/1000000, 2) for elem in strepto_gs]


florifile = os.path.join(p.parents[0], 'files', '06012020_floricoccus_database.tsv')
flori_gs = get_info(florifile, 'genome_length')
flori_gs = [round(elem/1000000, 2) for elem in flori_gs]

#%%

#Figure with all genome sizes combined (Floricoccus only rug plot)
fig = plt.figure(figsize=(8,6))
fig = sns.distplot(lac_gs, color = '#91cf60', label = 'Lactococcus', hist=0, kde=1, kde_kws = {'shade' : 1, 'linewidth' : 2})
fig = sns.distplot(strepto_gs, color = '#3288bd', label = 'Streptococcus', hist=0, kde=1, kde_kws = {'shade' : 1, 'linewidth' : 2})
fig = sns.distplot(flori_gs, color = '#b2182b', label = 'Floricoccus', hist=0, kde=0, rug =1, rug_kws={'linewidth' :2})

fig.legend(loc='upper right', frameon=0)
fig.set_title('Genome size distribution of the three analyzed genera', fontsize=16)
fig.set(xlabel='Genome size (in Mb)', ylabel='Density')
sns.set_style('whitegrid')
#remove dark grid from plot
# sns.set(style="whitegrid")
# fig.grid(False)
plt.tight_layout()
fig.figure.savefig(os.path.join(p.parents[0], 'figures', '270120_histogram_gs.png'), dpi=300, bbox_inches='tight')


#%%
######Make several subplots#######

#Make grid of plots (4 total) that can be filled
fig, axs = plt.subplots(2,2, figsize= (10,8)) #2x2 grid with determined figure size

#fill subplots with histograms
ax1 = sns.distplot(lac_gs, color = '#91cf60', label = 'Lactococcus', ax=axs[0,0], hist=0, kde=1, kde_kws = {'shade' : 1, 'linewidth' : 2})
ax1 = sns.distplot(strepto_gs, color = '#3288bd', label = 'Streptococcus', ax=axs[0,0], hist=0, kde=1, kde_kws = {'shade' : 1, 'linewidth' : 2})
ax1 = sns.distplot(flori_gs, color = '#b2182b', label = 'Floricoccus', ax=axs[0,0], hist=0, kde=0, rug =1, rug_kws={'linewidth' :2})
ax2 = sns.distplot(lac_gs, color = '#91cf60', label = 'Lactococcus', ax=axs[0,1], rug=1)
ax3 = sns.distplot(strepto_gs, color = '#3288bd', label = 'Streptococcus', ax=axs[1,0], rug=1)
ax4 = sns.countplot(flori_gs, color = '#b2182b', alpha= 0.5, label = 'Floricoccus', ax=axs[1,1])
#ax4 = sns.distplot(flori_gs,bins=20, color = '#b2182b', label = 'Floricoccus', ax=axs[1,1], kde=0, norm_hist=0, rug=1)
ax1.legend(loc='upper right')

#Titles
fig.suptitle('Genome sizes', fontsize=16, y=0.95)
ax1.set_title('Merged')
ax2.set_title('Lactococcus')
ax3.set_title('Streptococcus')
ax4.set_title('Floricoccus')

#Set x and y labels for all plots
# for ax in axs.flat:
    # ax.set(xlabel='Genome size (in Mb)', ylabel='Density')

ax1.set(xlabel='Genome size (in Mb)', ylabel='Density')
ax2.set(xlabel='Genome size (in Mb)', ylabel='Density')
ax3.set(xlabel='Genome size (in Mb)', ylabel='Density')
ax4.set(xlabel='Genome size (in Mb)', ylabel='Frequency')

#Sublabels for the plot
ax1.text(-0.05, 1.10, "A", ha = "left", va="top", transform=ax1.transAxes, size=12, weight = 'bold')
ax2.text(-0.05, 1.10, "B", ha = "left", va="top", transform=ax2.transAxes, weight = 'bold')
ax3.text(-0.05, 1.10, "C", ha = "left", va="top", transform=ax3.transAxes, weight = 'bold')
ax4.text(-0.05, 1.10, "D", ha = "left", va="top", transform=ax4.transAxes, weight = 'bold')

#Adjust layout to preserve title
plt.tight_layout()
plt.subplots_adjust(top=0.88)
fig.savefig(os.path.join(p.parents[0], 'figures', '230120_histogram_gs.png'), dpi=300)



