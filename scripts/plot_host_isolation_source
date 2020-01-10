#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 15:48:18 2020

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
florifile = os.path.join(p.parents[0], 'files', '06012020_floricoccus_database.tsv')
with open(florifile) as f:
    headers = f.readline().strip().split('\t')
    inds = {k: i for i, k in enumerate(headers)}
    hostname_i = inds['host_name']
    isolation_source_i = inds['isolation_source']
    
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

#Make dict containing all hosts as keys and values are dicts with all isolation sources as keys and 
#counts as values
host_isolations = {}

for k,v in isolations.items():
    count = {}
    for item in v:
        if item not in count:
            count[item] = 1
        else:
            count[item] += 1
    host_isolations[k] = count

df = pd.DataFrame(host_isolations)
df.T.plot(kind ='bar', legend = False, figsize = (15,5), width=2)

#plt.bar(range(len(host_isolations)), list(total_count.values()), align='center')
#plt.xticks(range(len(isolations)), isolations.keys(), rotation=90)
#plt.title("Genome sizes " + genus, fontsize = 14)
#plt.xlabel("Genome length\n(in Mb)", fontsize = 11)
#plt.ylabel("Number of species counted", fontsize = 11)
#plt.tight_layout()
host_count = {}       
for item in hosts:
    if item not in host_count and item != '':
        host_count[item] =1
    elif item !='':
        host_count[item] += 1
all_hosts = Counter(host_count).most_common()


top5 = Counter(host_count).most_common(5)
top5.sort()
others = [i for i in all_hosts if i not in top5]
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

#plt.bar(range(len(host_count)), numbers)
#plt.xticks(range(len(host_count)), labels, rotation=90)
#plt.title("Counted isolation sources related to hosts", fontsize= 14)
#plt.xlabel("Host")
#plt.ylabel("Amount of different isolation sources")

theme = plt.get_cmap('Blues')
colors = [theme(1. * i / len(numbers)) for i in range(len(numbers))]


fig1, ax1 = plt.subplots()
ax1.pie(numbers, labels=labels, autopct='%1.1f%%', colors=colors, startangle=90, radius = 1)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.title('Amount of isolation sources found in different hosts\n(Lactococcus)')
#plt.tight_layout()
#plt.savefig(os.path.join(p.parents[0], 'figures', '20200108_lactococcus_isolation_sources_per_host.tiff'), bbox_inches='tight')
