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
from brokenaxes import brokenaxes
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

def count(item_l):
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

#%% runcell 1
    
#host and associated sources in streptococcus 
#get info and the top5 occuring hosts with all isolation sources and a count of these

path = os.getcwd()
p = Path(path)

with open(os.path.join(p.parents[1], 'files', '20012020_streptococcus_database.tsv')) as f:
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
    host_isolations[host] = count(host_isolations[host])

#determine the top 5 most occuring hosts
count_hosts = [i for i in hosts if i != '']
counter_hosts = Counter(count_hosts)
top5_hosts = counter_hosts.most_common(5)    

#%% runcell 2 

#Plot the most occuring hosts and how often they appear in the dataset

h_labels= []
h_count = []
for i in top5_hosts:
    h_labels.append(i[0])
    h_count.append(i[1])

theme = plt.get_cmap('coolwarm')
colors = [theme(1. * i / len(h_labels)) for i in range(len(h_labels))]
 

#Make barplot    with broken y-axis    
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

# plot the same data on both axes
ax1.bar(range(len(h_count)), h_count, align='center', color= colors)
ax2.bar(range(len(h_count)), h_count, align='center', color= colors)

# zoom-in / limit the view to different portions of the data
ax1.set_ylim(7000, 7500)  # outliers only
ax2.set_ylim(0, 1000)  # most of the data

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
ax2.set_xticks(range(len(h_count)))
ax2.set_xticklabels(h_labels, rotation=90)
ax2.set_xlabel("Hosts")
ax2.set_ylabel("Frequency", y= 1.05)

ax1.set_title('The most common hosts of $Streptococcus$', fontsize=14, y=1.05) #$ around part that should be italics


#%% runcell 3 
   
#Plot top 10 sources of number one host

#get data (counts) and labels (isolation sources) in same order from dict count
human = host_isolations['Human, Homo sapiens']
for i in human:
    labels, data = keys_and_values(human)

#determine top 10 sources
most_common_10 = []
indexes = []
for i, count in enumerate(data):
    #fill top 10 with 10 source counts
    if len(most_common_10) < 10:
        most_common_10.append(count)
        indexes.append(i)
    #replace the minimum with the new count if its bigger
    else:
        min_count = min(most_common_10)
        ind_min = most_common_10.index(min_count)
        if count > min_count:
            most_common_10[ind_min] = count
            indexes[ind_min] = i

#Get labels
fig_labels=[]
for i in indexes:
    if labels[i] == '':
        fig_labels.append("Unknown")
    else:
        fig_labels.append(labels[i].capitalize())

theme = plt.get_cmap('Blues_r')
colors = [theme(1. * i / len(most_common_10)) for i in range(len(most_common_10))]
 

#Make barplot    with broken y-axis    
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

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

fig.savefig(os.path.join(p.parents[1], 'figures', '110220_isolation_sources_human.png'), dpi=300,bbox_inches='tight')

#%% runcell 4

'''
Try to use brokenaxes package, but problems with xticklabels....
'''

#make plot with broken y-axis
fig= plt.figure(figsize=(8,6))

bax= brokenaxes(ylims=((0,1100), (4000, 4600)), hspace=.15)


bax.bar(range(len(most_common_10)), most_common_10, align='center', color= colors)

bax.set_xticks(range(len(most_common_10)))
bax.set_xticklabels(fig_labels, rotation=90)

bax.set_xlabel('Isolation source')
bax.set_ylabel('Frequency')
fig.set_title('The most common isolation sources of $Streptococcus$ in Human', fontsize=16) #$ around part that should be italics

#remove dark grid from plot
# sns.set(style="whitegrid")
# fig.grid(False)
plt.tight_layout()
#fig.figure.savefig(os.path.join(p.parents[0], 'figures', '270120_histogram_gs.png'), dpi=300, bbox_inches='tight')   
   
  

           
          