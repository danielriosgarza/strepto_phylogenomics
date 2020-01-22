#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 10:53:23 2020

@author: meike
"""
# colours: spectral, coolwarm

import os
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter 
from scipy import stats


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

def histogram(item_l, title=None, xlabel=None, label=None, bins=10, color='#2166ac', rotation=0):
    '''
    Makes density histogram from given list.
    '''
    plt.hist(item_l, bins, density=1, alpha=0.75, color = color, label =label)
    plt.title(title, fontsize= 14)
    plt.xlabel(xlabel, fontsize= 12)
    plt.xticks(rotation=rotation)
    plt.ylabel("Density", fontsize= 12)
    plt.tight_layout() 
    
    #plt.savefig(os.path.join(p.parents[0], 'figures', )
    
    
    
path = os.getcwd()
p = Path(path)

#make histogram of genome sizes
lactofile = os.path.join(p.parents[0], 'files', '06012020_lactococcus_database.tsv')
lac_gs = get_info(lactofile, 'genome_length')
lac_gs = [round(elem/1000000, 2) for elem in lac_gs]
lac_hist = histogram(lac_gs,"Genome sizes of Lactococcus", "Genome sizes\n(in Mb)", bins=20, label = "Lactococcus")

lac_seqplatforms = get_info(lactofile, 'sequencing_platform')


streptofile = os.path.join(p.parents[0], 'files', '20012020_streptococcus_database.tsv')
strepto_gs = get_info(streptofile, 'genome_length')
strepto_gs = [round(elem/1000000, 2) for elem in strepto_gs]
strepto_hist = histogram(strepto_gs,"Genome sizes of Streptococcus", "Genome sizes\n(in Mb)", bins=20, color='#fee090', label = "Streptococcus")


florifile = os.path.join(p.parents[0], 'files', '06012020_floricoccus_database.tsv')
flori_gs = get_info(florifile, 'genome_length')
flori_gs = [round(elem/1000000, 2) for elem in flori_gs]

flori_hist = histogram(flori_gs,"Genome sizes of Floricoccus", "Genome sizes\n(in Mb)", bins=20, color='#b2182b', label= "Floricoccus")

plt.legend(loc='upper left', frameon=0)

############try to make breaks on y-axis
# f, (ax, ax2) = plt.subplots(2, 1, sharex=True)

# ax.plot(flori_hist)
# ax2.plot(flori_hist)

# # zoom-in / limit the view to different portions of the data
# ax.set_ylim(80, 100)  # outliers only
# ax2.set_ylim(0, 6)  # most of the data

# # hide the spines between ax and ax2
# ax.spines['bottom'].set_visible(False)
# ax2.spines['top'].set_visible(False)
# ax.xaxis.tick_top()
# ax.tick_params(labeltop=False)  # don't put tick labels at the top
# ax2.xaxis.tick_bottom()
#plt.ylim(0,10)
 

all_gs = {"Lactococcus" : lac_gs,
          "Streptococcus": strepto_gs,
          "Floricoccus" : flori_gs}
labels = []
data = []
for k,v in all_gs.items():
    data.append(k)
    labels.append(v)
    
histogram(all_gs, 'Genome sizes', 'Genomesizes\n(in Mb)', label = all_gs.keys())


plt.hist(all_gs)
   
    # platform_count ={}
    # for platform in seq_platforms:
    #     if platform not in platform_count:
    #         platform_count[platform] = 1
    #     else:
    #         platform_count[platform] += 1
        
    # top10 = Counter(platform_count).most_common(10)
    # top10.sort()
    
    # data = []
    # labels =[]
    # for item in top10:
    #     data.append(item[1])
    #     labels.append(item[0])




fig, ax = plt.subplots(figsize=(8,6))
ax.hist(lac_gs, 10, density=1, alpha=0.5, color = '#2166ac')
plt.title("Genome sizes of Lactococcus", fontsize= 14)
plt.xlabel("Genome sizes\n(in Mb)", fontsize= 12)
plt.ylabel("Density", fontsize= 12)
plt.tight_layout()
#plt.savefig(os.path.join(p.parents[0], 'figures', )

sns.distplot(lac_gs, kde=1, rug=1)
plt.xlabel("Genome size\n(in Mb)")
plt.xticks(rotation=45)

