#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 18:33:39 2020

@author: meike
"""
import numpy as np
import csv

#file on server: '/home/meiker/phylo_tree/roary/gff_files/orthology/gene_presence_absence.csv
file = '/home/meike/Downloads/gene_presence_absence.csv'

#Look at presence of genes (save as binary)
genes2bin = {}
gene2ids = {}
with open(file) as csv_file:
    f = csv.reader(csv_file, delimiter=',')
    for ind, line in enumerate(f):
        if ind == 0:
            ids = line[14::]
            pass
        else:
            gene2ids[line[0]] = []
            binary = np.zeros(len(ids))
            for i, pres in enumerate(line[15::]):
                if pres != '':
                    binary[i] = 1
                    gene2ids[line[0]].append(ids[i])
            genes2bin[line[0]] = binary 
        
#look for core and extended core genes (all above 0.95 presence)            
core = []

for k in genes2bin:
    score = np.count_nonzero(genes2bin[k])/len(genes2bin[k])
    if round(score,2) >= 0.95:
        core.append(k)

missing = {}
for g in core:
    if len(gene2ids[g]) != 204:
        for item in ids:
            if item not in gene2ids[g]:
                if item not in missing:
                    missing[item] = 1
                else:
                    missing[item] += 1