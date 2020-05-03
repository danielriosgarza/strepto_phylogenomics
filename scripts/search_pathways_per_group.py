#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:19:54 2020

@author: meike
"""

'''
Look for pathways with more or less genes present per group. Looks for genes with p-value <0.05 and sets names of genes in a list. These can be used on the website Enrichr to look for affected pathways.
'''

import os
from pathlib import Path
from datetime import date
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import scipy.stats as sts
import seaborn as sns


path = os.getcwd()
p = Path(path)

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)

bin_greater = os.path.join(p.parents[0], 'files', 'binary_table', '24042020_pvalues_genes_greater.tsv')

bin_less = os.path.join(p.parents[0], 'files', 'binary_table', '24042020_pvalues_genes_less.tsv')

group2grless2genes = {i : { 'greater': [], 'less' : []} for i in range(8)}

#determine number of cols and save the lines
with open(bin_greater) as f:
    f.readline()
    for line in f:
        a = line.strip().split('\t')
        gene_num = a[0].split('_')[1]
        for ind, value in enumerate(a[1::]):
            if float(value) <= 0.05:
                group2grless2genes[ind]['greater'].append(gene_num)
                
with open(bin_less) as f:
    f.readline()
    for line in f:
        a = line.strip().split('\t')
        gene_num = a[0].split('_')[1]
        for ind, value in enumerate(a[1::]):
            if float(value) <= 0.05:
                group2grless2genes[ind]['less'].append(gene_num)
                
#Dict: group --> less/greater --> gene names, write lists into files to feed to Enrichr
g_names_needed = []
for gr in group2grless2genes:
    for di in group2grless2genes[gr]:
        for num in group2grless2genes[gr][di]:
            if num not in g_names_needed:
                g_names_needed.append(num)
g_names_needed = sorted(g_names_needed)

genes = {}
with open(os.path.join(p.parents[0], 'files', 'binary_table', '14042020_binary_table_sorted.tsv')) as f:
    f.readline()
    lines = [l for l in f]

for ind in g_names_needed:
    ge = lines[int(ind)].strip().split('\t')[2]
    ge = ge.split(',')
    genes[ind] = ge
        

groups_genes = {}   

for gr in group2grless2genes:
    groups_genes[gr] = {}
    for Dir in group2grless2genes[gr]:
        groups_genes[gr][Dir] = []
        for num in group2grless2genes[gr][Dir]:
            name = genes[num]
            if 'N/A' not in name:
                groups_genes[gr][Dir] += name
                
# for i_gr, gr in enumerate(groups_genes):
#     folder = os.path.join(p.parents[0], 'files', 'pathway_search','group' + str(i_gr + 1))
#     if not os.path.isdir(folder):
#         os.mkdir(folder)
#     for grless in groups_genes[gr]:
#         with open(folder + '/' + grless + '_gene_list.txt', 'w') as f:
#             for ge_name in groups_genes[gr][grless]:
#                 f.write(ge_name + '\n')


