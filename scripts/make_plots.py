#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 12:24:35 2020

@author: meike
"""

'''
Make Plots
'''
#imports and funcitons

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.gridspec as gridspec
from collections import Counter
import heapq

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
                        item_l.append(str(item))
    return item_l

def get_count(item_l):
    '''
    Counts items in a list, returns lists containg data and labels.
    '''    
    count = {}
    for item in item_l:
        if item not in count:
            count[item] = 1
        else:
            count[item] += 1
    return count
    
def keys_and_values(dict_):
    '''
    Needs a dict and returns keys and values as lists
    '''    
    keys = []
    values = []
    for k,v in dict_.items():
        keys.append(k)
        values.append(v)
    return keys, values

def determine_topN(countD, Nitems=10):
    '''
    Needs a count dict with keys labels and values the according count.
    Uses the function keys_and_values!
    Returns to lists (labels and counts) of the top N counted items of the dict.
    Default of N is 10.
    '''
    #get data (counts) and labels in same order from dict count
    labels, data = keys_and_values(countD)
    
    #determine top N sources
    topN = []
    indexes = []
    for i, count in enumerate(data):
        #fill topN with N source counts
        if len(topN) < Nitems:
            topN.append(count)
            indexes.append(i)
        #replace the minimum with the new count if its bigger
        else:
            min_count = min(topN)
            ind_min = topN.index(min_count)
            if count > min_count:
                topN[ind_min] = count
                indexes[ind_min] = i
                
    #Get labels
    fig_labels=[]
    for i in indexes:
        if labels[i] == '':
            fig_labels.append("Unknown")
        else:
            fig_labels.append(labels[i].capitalize())
    
    return fig_labels, topN

def simple_barplot(data, labels, color='coolwarm'):
    '''
    Simple barplot w/o axis labels or title. Neeeds two lists with numerical data and labels (for xticks).
    '''
    theme = plt.get_cmap(color)
    plabels, pdata = determine_topN(pig_sources, 5)
    colors = [theme(1. * i / len(data)) for i in range(len(data))]

    fig, ax = plt.subplots()
    
    ax.bar(range(len(data)), data, align='center', color= colors)
    ax.set_xticks(range(len(data)))
    ax.set_xticklabels(labels, rotation=90)
    return fig, ax

def simple_pie(data, labels, color='coolwarm', origin =0):
    '''
    '''
    fig, ax = plt.subplots()
    #autopct='%1.1f%%' gives percentages to pie, pctdistance changes postition of percentages
    ax.pie(data, autopct='%1.1f%%', colors=colors, startangle=origin, radius = 1, pctdistance=1.15) 
    
    #loc in combi w/ bbocx --> loc tells matplotlib which part of bounding box should be placed at the arguments of bbox_to_anchor
    ax.legend(labels, loc= "center left", bbox_to_anchor=(0.9, 0.5))
    
    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    return fig, ax
    
#%% runcell 1
    
#host and associated sources in streptococcus 
#get info and the top5 occuring hosts with all isolation sources and a count of these

path = os.getcwd()
p = Path(path)

with open(os.path.join(p.parents[0], 'files', '20012020_streptococcus_database.tsv')) as f:
          headers = f.readline().strip().split('\t')
          inds = {k: i for i, k in enumerate(headers)}
          hostname_i = inds['host_name']
          isolation_source_i = inds['isolation_source']
          #make lists of hosts and the related isolation source
          hosts = []
          source=[]
          for line in f:
              a = line.strip().split('\t')
              hosts.append(a[hostname_i]) 
              source.append(a[isolation_source_i]) 

host_isolations = {k :[] for k in hosts if k!=''}

for i, item in enumerate(hosts):
    if item != '':
        host_isolations[item].append(source[i])

#Make nested dict with hosts as keys and count dict of the sources       
for host in host_isolations:
    host_isolations[host] = get_count(host_isolations[host])

#determine the top 5 most occuring hosts
count_hosts = [i for i in hosts if i != '']
counter_hosts = Counter(count_hosts)
top5_hosts = counter_hosts.most_common(5)    


#%% runcell 2 

#Plot the most occuring hosts and how often they appear in the dataset (barplot)

h_labels= []
h_count = []
for i in top5_hosts:
    h_labels.append(i[0])
    h_count.append(i[1])

theme = plt.get_cmap('coolwarm')
colors = [theme(1. * i / len(h_labels)) for i in range(len(h_labels))]
 

#Make barplot with broken y-axis    
fig1, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

#plot the same data on both axes
ax1.bar(range(len(h_count)), h_count, align='center', color= colors)
ax2.bar(range(len(h_count)), h_count, align='center', color= colors)

#zoom-in / limit the view to different portions of the data
ax1.set_ylim(7000, 7500)  # outliers only
ax2.set_ylim(0, 1000)  # most of the data

#hide the spines between ax1 and ax2
ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)

#ax1.xaxis.tick_top()
ax1.tick_params(labeltop=False)  #don't put tick labels at the top
ax2.xaxis.tick_bottom() #make ticks at the bottom of x-axis

#cut-out diagonal lines
d = .015  #how big to make the diagonal lines in axes coordinates
kwargs = dict(transform=ax1.transAxes, color='lightgray', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        #top-left diagonal
kwargs.update(transform=ax2.transAxes)  #switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  #bottom-left diagonal

#make labels
ax2.set_xticks(range(len(h_count)))
ax2.set_xticklabels(h_labels, rotation=90)
ax2.set_xlabel("Hosts")
ax2.set_ylabel("Frequency", y= 1.05)

ax1.set_title('The most common hosts of $Streptococcus$', fontsize=14, y=1.05) #$ around part that should be italics

#%% runcell 3 

#Pieplot: Plot the most occuring hosts and how often they appear in the dataset


fig2, ax = plt.subplots()
#autopct='%1.1f%%' gives percentages to pie, pctdistance changes postition of percentages
ax.pie(h_count, autopct='%1.1f%%', colors=colors, startangle=0, radius = 1, pctdistance=1.15) 

#loc in combi w/ bbocx --> loc tells matplotlib which part of bounding box should be placed at the arguments of bbox_to_anchor
ax.legend(h_labels, loc= "center left", bbox_to_anchor=(0.9, 0.5))

ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
#pad defines location of title to plot
plt.title('The most common hosts of $Streptococcus$', pad=1)

#plt.savefig(os.path.join(savedir), bbox_inches='tight')
    


#%% runcell 4
   
#Plot top 10 sources of number one host

#get data (counts) and labels (isolation sources) in same order from dict count
fig_labels, most_common_10 = determine_topN(host_isolations['Human, Homo sapiens'])

theme = plt.get_cmap('coolwarm')
colors = [theme(1. * i / len(most_common_10)) for i in range(len(most_common_10))]
 

#Make barplot    with broken y-axis    
fig3, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

# plot the same data on both axes
ax1.bar(range(len(most_common_10)), most_common_10, align='center', color= colors)
ax2.bar(range(len(most_common_10)), most_common_10, align='center', color= colors)

# zoom-in / limit the view to different portions of the data
ax1.set_ylim(4000, 4600)  # outliers only
ax2.set_ylim(0, 1100)  # most of the data

# hide the spines between ax1 and ax2
ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)

#ax1.xaxis.tick_top()
ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

#cut-out diagonal lines
d = .015  #how big to make the diagonal lines in axes coordinates
kwargs = dict(transform=ax1.transAxes, color='lightgray', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal

#make labels
ax2.set_xticks(range(len(most_common_10)))
ax2.set_xticklabels(fig_labels, rotation=90)
ax2.set_xlabel("Isolation sources")
ax2.set_ylabel("Frequency", y= 1.05)

ax1.set_title('The most common isolation sources of $Streptococcus$ in Human', fontsize=14, y=1.05) #$ around part that should be italics

#fig.savefig(os.path.join(p.parents[1], 'figures', '110220_isolation_sources_human.png'), dpi=300,bbox_inches='tight')

#%% runcell 5

#Look at other hosts
#pig
#combine all Tonsil sources
pig_sources ={}
for k,v in host_isolations['Pig, Sus scrofa'].items():
    if "onsil" in k:
        if "Tonsil" not in pig_sources:
            pig_sources['Tonsil'] = v
        else:
            pig_sources['Tonsil']+=v
    else:
        pig_sources[k] = v
    
plabels, pdata = determine_topN(pig_sources, 5)
colors = [theme(1. * i / len(pdata)) for i in range(len(pdata))]

fig4, ax = plt.subplots()

ax.bar(range(len(pdata)), pdata, align='center', color= colors)
ax.set_xticks(range(len(pdata)))
ax.set_xticklabels(plabels, rotation=90)
ax.set_xlabel("Isolation sources")
ax.set_ylabel("Frequency")
ax.set_title('The most common isolation sources of $Streptococcus$ in Pig', fontsize=14, y=1.05)

#Cow
cow_sources ={}
for k,v in host_isolations['Cow, Bos taurus'].items():
    if "ilk" in k:
        if "Milk" not in cow_sources:
            cow_sources['Milk'] = v
        else:
            cow_sources['Milk']+=v
    elif "available" in k:
        k=''
        if '' not in cow_sources:
            cow_sources[''] = v
        else:
            cow_sources[''] += v
    else:
        cow_sources[k] = v

clabels, cdata = determine_topN(cow_sources, 5)
#fig4, ax = simple_pie(cdata, clabels, origin=90)
fig5, ax = simple_barplot(cdata, clabels)
ax.set_xlabel("Isolation sources")
ax.set_ylabel("Frequency")
ax.set_title('The most common isolation sources\n of $Streptococcus$ in Cow', fontsize=14, y=1.05)

#Horse
horse_sources ={}
for k,v in host_isolations['Horse, Equus caballus'].items():
    if "Nasal" in k or "Naso" in k:
        if "Nasopharynx" not in horse_sources:
            horse_sources['Nasopharynx'] = v
        else:
            horse_sources['Nasopharynx']+=v
    elif 'bscess' in k:
        if "Abscess" not in horse_sources:
            horse_sources['Abscess'] = v
        else:
            horse_sources['Abscess']+=v
    else:
        horse_sources[k] = v

holabels, hodata = determine_topN(horse_sources, 5)
#fig4, ax = simple_pie(cdata, clabels, origin=90)
fig6, ax = simple_barplot(hodata, holabels)
ax.set_xlabel("Isolation sources")
ax.set_ylabel("Frequency")
ax.set_title('The most common isolation sources\n of $Streptococcus$ in Horse', fontsize=14, y=1.05)


#Fish Oreochromis niloticus
fish_sources ={}
for k,v in host_isolations['Oreochromis niloticus'].items():
    if "disease" in k or 'infect' in k:
        if "Diseased Fish" not in fish_sources:
            fish_sources["Diseased Fish"] = v
        else:
            fish_sources["Diseased Fish"] += v
    else:
        fish_sources[k] = v

flabels, fdata = determine_topN(fish_sources, 5)
#fig4, ax = simple_pie(cdata, clabels, origin=90)
fig7, ax = simple_barplot(fdata, flabels)
ax.set_xlabel("Isolation sources")
ax.set_ylabel("Frequency")
ax.set_title('The most common isolation sources\n of $Streptococcus$ in Fish', fontsize=14, y=1.05)


#%% runcell 6

#Combine plots in a single figure

figall, axs = plt.subplots(3,2, figsize= (10,8)) #2x2 grid with determined figure size

#fill subplots 
ax1 = axs[0,0]
axs[0,1] 
axs[1,0] 
axs[1,1] 
axs[2,0] 
axs[2,1] 

#Titles
#figall.suptitle('Genome sizes', fontsize=16, y=0.95)
# ax1.set_title('Merged')
# ax2.set_title('Lactococcus')
# ax3.set_title('Streptococcus')
# ax4.set_title('Floricoccus')

#Set x and y labels for all plots
# for ax in axs.flat:
    # ax.set(xlabel='Genome size (in Mb)', ylabel='Density')

#Sublabels for the plot
#axs[0,0].text(-0.05, 1.10, "A", ha = "left", va="top", size=12, weight = 'bold')
# axs[0,1].text(-0.05, 1.10, "B", ha = "left", va="top", transform=axs[0,1].transAxes, weight = 'bold')
# axs[1,0].text(-0.05, 1.10, "C", ha = "left", va="top", transform=axs[1,0].transAxes, weight = 'bold')
# axs[1,1].text(-0.05, 1.10, "D", ha = "left", va="top", transform=axs[1,1].transAxes, weight = 'bold')
# axs[2,0].text(-0.05, 1.10, "E", ha = "left", va="top", transform=axs[2,0].transAxes, weight = 'bold')
# axs[2,1].text(-0.05, 1.10, "F", ha = "left", va="top", transform=axs[2,1].transAxes, weight = 'bold')

#Adjust layout to preserve title
plt.tight_layout()
plt.subplots_adjust(top=0.88)
#fig.savefig(os.path.join(p.parents[0], 'figures', '230120_histogram_gs.png'), dpi=300)   
  

           
          