#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:42:19 2020

@author: meike
"""


'''
Binary table sorting.

The binary table will be sorted according to the appearence of the orthologues (appearence of 1).

Pan-genome : presence of orthologue in genomes
Core : 100%
Extended core : 90-100%
Shell : 15-90%
Cloud : <15%
'''

import os
from pathlib import Path
from datetime import date
import numpy as np
    
path = os.getcwd()
p = Path(path)

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)

bin_file = os.path.join(p.parents[0], 'files', 'binary_table', '01042020_binary_table_prep.tsv')

#determine number of cols and save the lines
with open(bin_file) as f:
    lines = [line for line in f]

ncols = len(lines[0].split('\t'))
    

#set the binary part into numpy array
data = np.loadtxt(bin_file, delimiter = '\t', skiprows = 1, usecols = range(3, ncols))
    
#determine the number of ones per line (appearence of gene in the genomes)
scores = np.count_nonzero(data, axis=1)

#determine the percentage of the gene presence
pscores = scores/(ncols - 1)

#sort the scores in descending (gives indexes of the scores)
sorted_scores = np.argsort(pscores)[::-1]

    
with open(os.path.join(p.parents[0], 'files', 'binary_table', today + '_binary_table_sorted.tsv'), 'w') as f:
    f.write('Pan-genome\tAppearance (in %)\t' + lines[0])
    for i in sorted_scores:
        #because the header is missing index of line is +1
        i += 1 
        perc = round(pscores[i]*100,2)
        if pscores[i] == 1:
            f.write('Core\t' + str(perc) + '\t' + lines[i])
        elif pscores[i] >= 0.9:
            f.write('Extended Core\t' + str(perc) + '\t' + lines[i])
        elif 0.9 > pscores[i] >= 0.15:
            f.write('Shell\t' + str(perc) + '\t' + lines[i])
        elif 0.15 > pscores[i]:
            f.write('Cloud\t' + str(perc) + '\t' + lines[i])
            

  