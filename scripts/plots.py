#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 16:31:04 2019

@author: meike
"""
'''Tests to visualize data'''

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

data = pd.read_csv('/home/meike/strepto_phylogenomics/files/strepto_genomes_quality.tsv', delimiter="\t") 
print(data.head(2))
test = data.head(20)

#data.plot(x ='genome.genome_name', y='genome.genome_length', kind='hist')
#data.plot(x ='genome.genome_name', y='genome.host_name', kind='box')

#sns.countplot(y='genome.host_name', hue='genome.genome_name', data=data, palette="Greens_d")

#df=test.groupby(['genome.genome_name','genome.host_name']).size()
#df=df.unstack()
#df.plot(kind='bar')

countries = {}

for item in data['genome.isolation_country']:
    if item not in countries:
        countries[item] = 1
    else:
        countries[item] += 1
sources = {}        
for i in data['genome.isolation_source']:
    if i not in sources: 
        sources[i] = 1
    else:
        sources[i] += 1
top_10_sources = dict(Counter(sources).most_common(10))
plt.bar(range(len(top_10_sources)), top_10_sources.values())
plt.xticks(range(len(top_10_sources)), top_10_sources.keys())
plt.show()

a = dict(Counter(countries).most_common(10))
print(a)

plt.bar(range(len(a)), list(a.values()))
plt.xticks(range(len(a)), list(a.keys()))
plt.show()

hosts= []
typus =[]
for col in test[['genome.host_name']]:
    hosts.append(test[col])
for j in test[['genome.oxygen_requirement']]:
    typus.append(test[j])
combi =[]
for a in range(len(hosts)):
    print(a)
    #combi.append(a, b)
    
print(combi)
#plt.plot(hosts) and plt.plot(typus)
#plt.bar(range(len(hosts)), hosts)
#plt.show

