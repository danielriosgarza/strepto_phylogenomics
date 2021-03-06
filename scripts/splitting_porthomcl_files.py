#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 16:32:22 2020

@author: meike
"""

'''
Split Ortholouge/Paralogue search files
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


original_obash = os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', '24022020_orthologs.sh')
o_bash_dir = os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl')

p_bash = os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', '24022020_paralogs.sh')

indexes = []
#look how many orthologues were already found and remove them from the file before splitting it
with open(original_obash) as f:
    last = False
    i = 1
    for line in f:
        i += 1
        if i == 314: #change number to amount of files in orthology dir
            last = True
        if last:
            line = line.strip().split(' ')
            indexes.append(line[-1])

size = round(len(indexes)/10)

line = 'porthomclPairsOrthologs.py -t /home/meiker/orthomcl/taxon_list -b /home/meiker/orthomcl/besthit -o /home/meiker/orthomcl/orthologs -x '
for i in range(1, 11):
    with open(o_bash_dir + '/orthologs/' + today + '_orthologs' + str(i) + '.sh', 'w') as f:  
        print(i)
        if i == 1:
            numbers = indexes[:size]
            for n in numbers:
                f.write(line + n + '\n')        
        if i == 2:
            numbers = indexes[size:size*2]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 3:
            numbers = indexes[size*2:size*3]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 4:
            numbers = indexes[size*3:size*4]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 5:
            numbers = indexes[size*4:size*5]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 6:
            numbers = indexes[size*5:size*6]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 7:
            numbers = indexes[size*6:size*7]
            for n in numbers:
                f.write(line + n + '\n')   
        if i == 8:
            numbers = indexes[size*7:size*8]
            for n in numbers:
                f.write(line + n + '\n')   
        if i == 9:
            numbers = indexes[size*8:size*9]
            for n in numbers:
                f.write(line + n + '\n')   
        if i == 10:
            numbers = indexes[size*9::]
            for n in numbers:
                f.write(line + n + '\n')   


p_bash = os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', '24022020_paralogs.sh')

indexes = []
#get the indexes to split them over new files
with open(p_bash) as f:
    for line in f:
        line = line.strip().split(' ')
        indexes.append(line[-1])

size = round(len(indexes)/12)

line = 'porthomclPairsInParalogs.py -t /home/meiker/orthomcl/taxon_list -q /home/meiker/orthomcl/paralogTemp -o /home/meiker/orthomcl/ogenes -p /home/meiker/orthomcl/paralogs -x '

for i in range(1, 13):
    with open(os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'paralogs', today + '_paralogs' + str(i) + '.sh'), 'w') as f:
        if i == 1:
            numbers = indexes[:size]
            for n in numbers:
                f.write(line + n + '\n')        
        if i == 2:
            numbers = indexes[size:size*2]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 3:
            numbers = indexes[size*2:size*3]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 4:
            numbers = indexes[size*3:size*4]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 5:
            numbers = indexes[size*4:size*5]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 6:
            numbers = indexes[size*5:size*6]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 7:
            numbers = indexes[size*6:size*7]
            for n in numbers:
                f.write(line + n + '\n')   
        if i == 8:
            numbers = indexes[size*7:size*8]
            for n in numbers:
                f.write(line + n + '\n')   
        if i == 9:
            numbers = indexes[size*8:size*9]
            for n in numbers:
                f.write(line + n + '\n')   
        if i == 10:
            numbers = indexes[size*9:size*10]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 11:
            numbers = indexes[size*10:size*11]
            for n in numbers:
                f.write(line + n + '\n')
        if i == 12:
            numbers = indexes[size*11::]
            for n in numbers:
                f.write(line + n + '\n')
     
        