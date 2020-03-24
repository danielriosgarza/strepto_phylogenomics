#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 08:29:40 2020

@author: meike
"""

import os
from pathlib import Path
from datetime import date


    
path = os.getcwd()
p = Path(path)

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)


original_obash = os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', '24022020_orthologs.sh')

indexes = []
with open (original_obash) as f:
    for line in f:
        line = line.strip().split(' ')
        indexes.append(line[-1])

index2id = {}

with open(os.path.join(p.parents[0], 'files', 'taxon_list')) as f:
    i = 1
    for line in f:
        line = line.strip()
        index2id[str(i)] = line
        i += 1

done = os.listdir('/home/meiker/orthomcl/orthologs')

ids = [] 
for i in done:
    i = i.split('.')[0]
    ids.append(i)

indexes_missing = []
for i in indexes:
    _id = index2id[i]
    if _id not in ids:
        indexes_missing.append(i)

size = round(len(indexes_missing)/10)

line = 'porthomclPairsOrthologs.py -t /home/meiker/orthomcl/taxon_list -b /home/meiker/orthomcl/besthit -o /home/meiker/orthomcl/orthologs -x '
for i in range(1,11):        
    with open (os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'orthologs') + '/' + today + '_orthologs_missing' + str(i) + '.sh', 'w') as f:
        if i == 1:
            numbers = indexes_missing[:size]
            for n in numbers:
                f.write(line + n + '\n')        
        if i == 2:
            numbers = indexes_missing[size:size*2]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 3:
            numbers = indexes_missing[size*2:size*3]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 4:
            numbers = indexes_missing[size*3:size*4]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 5:
            numbers = indexes_missing[size*4:size*5]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 6:
            numbers = indexes_missing[size*5:size*6]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 7:
            numbers = indexes_missing[size*6:size*7]
            for n in numbers:
                f.write(line + n + '\n')   
        if i == 8:
            numbers = indexes_missing[size*7:size*8]
            for n in numbers:
                f.write(line + n + '\n')   
        if i == 9:
            numbers = indexes_missing[size*8:size*9]
            for n in numbers:
                f.write(line + n + '\n')   
        if i == 10:
            numbers = indexes_missing[size*9::]
            for n in numbers:
                f.write(line + n + '\n') 
