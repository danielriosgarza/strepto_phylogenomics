#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 10:49:46 2020

@author: meike
"""

'''
Check the 1Mb genomes for core genes in order to investigate completness of the sequence (incomplete genome or genome reduction?)
'''

import os
from pathlib import Path

path = os.getcwd()
p = Path(path)

bin_file = os.path.join(p.parents[0], 'files', 'binary_table', '14042020_binary_table_sorted.tsv')

#save for each strain the number of (extended) core genes it has to determine what the average is
core = {} 
num_genes = 0
with open(bin_file) as f:
    ids = f.readline().strip().split('\t')[5::]
    for i in ids:
        core[i] = 0
    for line in f:
        a = line.strip().split('\t')
        if a[0] == 'Core' or a[0] == 'Extended Core':
            num_genes += 1
            for j, num in enumerate(a[5::]):
                if num == '1':
                    core[ids[j]] += 1

#calculate how many of the core genes each strain has
counts = [] 

for k,v in core.items():
    counts.append(v/num_genes)
    
average = sum(counts)/len(counts)
print('Average: ', average)

#check the number of core genes of the smaller genomes (ids from tree)
check_ids = ['streptococcus_11898','streptococcus_11946','streptococcus_11955','streptococcus_11902','streptococcus_11956']

indexes = [ids.index(ind) for ind in check_ids]
for index, i in enumerate(indexes):
    print('ID: ', check_ids[index], 'Core genes: ', counts[i])
    
    