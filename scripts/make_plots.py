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
import seaborn as sns
import pandas as pd
import matplotlib.gridspec as gridspec
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
    
    


#%%
#host and associated sources

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
        

           
          