#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:42:19 2020

@author: meike
"""


'''
Binary table sorting.

The binary table will be sorted according to the appearence of the orthologues (appearence of 1).

Pan-genome part : presence of orthologue in genomes
Core : 100%
Extended core : 90-100%
Shell : 15-90%
Cloud : <15%
'''

import os
from pathlib import Path
from datetime import date
    
path = os.getcwd()
p = Path(path)

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)

    
core = []
extended_core = []
shell = []
cloud = []

#sorted(input().split(),key=lambda x:-bin(int(x)).count("1")))
all_scores = []
    
with open(os.path.join(p.parents[1], 'tests', 'Files', 'test_binary.tsv')) as f:
        f.readline()
        for i, line in enumerate(f):
            i += 1
            numbers = line.strip().split('\t')[1::]
            score = numbers.count('1')/len(numbers)
            all_scores.append((score, i))
            
sorted_scores = sorted(all_scores, key=lambda x: x[0])    

#file into np array and look for indexes to write in correct order into new file      