#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 16:10:12 2020

@author: meike
"""

import matplotlib.pyplot as plt
import os
from pathlib import Path
import pandas as pd
import numpy as np
from collections import Counter 

path = os.getcwd()
p = Path(path)


lactofile = os.path.join(p.parents[0], 'files', '06012020_lactococcus_database.tsv')

def find_iso_sources_hosts(infile, genus):
    '''
    '''
    with open(infile) as f:
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
    
    isolations = {k : [] for k in hosts if k != ''}
    
    for i in range(len(hosts)):
        if hosts[i] != '':
            isolations[hosts[i]].append(source[i])
    
    #count the different isolations sources that are associated with certain host
    host_count = {}       
    for k,v in isolations.items():
        host_count[k] = len(set(v))
    
    top5 = Counter(host_count).most_common(5)
    top5.sort()
    
########Look at it again!!!!!!!!!!!!!!!! not sure f it works ...confused ....AHHHHHH!!!! ####################
    #ohh and check and write rest of the function
    others = [i for i in host_count if set(i) not in top5]
    rest = {}
    for i in range(len(others)):
        if others[i][1] not in rest:
            rest[others[i][1]] = 1
        else:
            rest[others[i][1]] += 1
    
    top5.append(('other', len(rest)))   
        
    top5.sort()    
        
    numbers = []
    labels =[]
    for item in top5:
        numbers.append(item[1])
        labels.append(item[0])
    
    theme = plt.get_cmap('Blues')
    colors = [theme(1. * i / len(numbers)) for i in range(len(numbers))]
    
    
    fig1, ax1 = plt.subplots()
    ax1.pie(numbers, labels=labels, autopct='%1.1f%%', colors=colors, startangle=90, radius = 1)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.title('Amount of isolation sources found in different hosts\n(Lactococcus)')

    #plt.savefig(os.path.join(p.parents[0], 'figures', '20200108_lactococcus_isolation_sources_per_host.tiff'), bbox_inches='tight')
