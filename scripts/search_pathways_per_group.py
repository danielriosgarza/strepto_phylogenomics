#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:19:54 2020

@author: meike
"""


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
                
# with open('/home/meike/tests/Files/interesting_genes.tsv', 'w') as f:
#     f.write ('greater\tless\n')
#     for i, v in enumerate(group2grless2genes[0]['greater']):
#         f.write(v + '\t' + group2grless2genes[0]['less'][i] + '\n')
            


