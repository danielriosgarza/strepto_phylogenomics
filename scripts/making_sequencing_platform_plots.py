#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 13:46:34 2020

@author: meike
"""


#fig, ax = plt.subplots(4, 4) --> makes 4x4 subplots in one pic

import matplotlib.pyplot as plt
import os
from pathlib import Path
from collections import Counter 


def seq_platform_plot(file, genus):
    '''
    Makes bar plot of top 10 used sequening platforms. x_ax = platforms, y_ax=times counted
    '''
    seq_platforms = []
    with open (file) as f:
        headers = f.readline().strip().split('\t')
        headers_inds = {name: i for i, name in enumerate(headers)}
        
        for line in f:
            a = line.strip().split('\t')
            for i, item in enumerate(a):
                if headers_inds['sequencing_platform'] == i and item != '':
                    seq_platforms.append(item)
    
    platform_count ={}
    for platform in seq_platforms:
        if platform not in platform_count:
            platform_count[platform] = 1
        else:
            platform_count[platform] += 1
        
    top10 = Counter(platform_count).most_common(10)
    top10.sort()
    
    data = []
    labels =[]
    for item in top10:
        data.append(item[1])
        labels.append(item[0])
    
        
    plt.bar(range(len(labels)), list(data), align='center')
    plt.xticks(range(len(labels)), labels, rotation=45, ha ='right')
    plt.title("Top 10 Sequencing Platforms used for " + genus, fontsize = 14)
    plt.xlabel("Sequencing Platform", fontsize = 11)
    plt.ylabel("Times mentioned", fontsize = 11)
    plt.tight_layout()
    # plt.savefig(saving_dir)

path = os.getcwd()
p = Path(path)

lactofile = os.path.join(p.parents[0], 'files', '06012020_lactococcus_database.tsv')
seq_platform_plot(lactofile, "Lactococcus")