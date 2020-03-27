#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 11:55:25 2020

@author: meike
"""


'''
Prepare Binary table of Pan-genome. Filter out the paralogs from the orthologs file.
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


#make dict containing all paralogs that they can be removed from the orthologs file
paralogs = {}
with open(os.path.join(p.parents[0], 'files', 'binary_table', 'all.par.group')) as f:
    for line in f:
        a = line.strip().split('\t')
        for i in a:
            _id = i.split('|')[0]
            para = i.split('|')[1]
            if _id not in paralogs:
                paralogs[_id] = [para]
            else:
                paralogs[_id] += [para]

i = 0 
      
with open (os.path.join(p.parents[0], 'files', 'binary_table', 'all.ort.group')) as f:
    with open(os.path.join(p.parents[0], 'files', 'binary_table', today + '_one1one_all.ort.group'), 'w') as f2:
        first = f.readline()
        for line in f:
            a = line.strip().split('\t')
            double = []
            ids = []
            orthos = []
            orthogroup = ''
            for pair in a:
                id_ = pair.split('|')[0]
                ortho = pair.split('|')[1]
                if id_ not in ids:
                    ids.append(id_)
                    orthos.append(ortho)
                else:
                    double.append(id_)
            for i, _id in enumerate(ids):
                if _id not in double:
                    orthogroup += _id + '|' + orthos[i] + '\t'
            orthogroup += '\n'
            f2.write(orthogroup)
                
            