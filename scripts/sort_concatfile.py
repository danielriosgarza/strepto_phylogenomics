#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 08:47:26 2020

@author: meike
"""


'''
Sort the concat file: genomes with known species first followed by unclassified --> ensures that these are taken into the tree, when redundant seqs are removed (removing_similar_seqs_from_concat_alignment.py)
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


concatfile = '/home/meiker/phylo_tree/iqtree/concat_alignments'
#concatfile =  '/home/meike/strepto_phylogenomics/files/random_files/concatenated_test'

seqs = {}
ids = []
with open(concatfile) as f:
    name = ''
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            name = line[1:]
            seqs[name] = ''
            ids.append(name)
        else:
            seqs[name] += line
            
id2species = {}

with open(os.path.join(p.parents[0], 'files', '03032020_streptococcus_database_final.tsv')) as f:
    headers = f.readline().strip().split('\t')
    species_ind = headers.index('species')
    for line in f:
        a = line.strip().split('\t')
        species = a[species_ind]
        id_ = a[0]
        id2species[id_] = species

with open(os.path.join(p.parents[0], 'files', '06012020_floricoccus_database.tsv')) as f:
    headers = f.readline().strip().split('\t')
    species_ind = headers.index('species')
    for line in f:
        a = line.strip().split('\t')
        species = a[species_ind]
        id_ = a[0]
        id2species[id_] = species

with open(os.path.join(p.parents[0], 'files', '06012020_lactococcus_database.tsv')) as f:
    headers = f.readline().strip().split('\t')
    species_ind = headers.index('species')
    for line in f:
        a = line.strip().split('\t')
        species = a[species_ind]
        id_ = a[0]
        id2species[id_] = species

sorted_ids = []
unknown = []
for id_ in ids:
    species = id2species[id_]
    if any(s.isdigit() for s in species) or 'uncultured' in species:
        unknown.append(id_)
    else:
        sorted_ids.append(id_)

for num in unknown:
    sorted_ids.append(num)
    
with open('/home/meiker/phylo_tree/roary/' + today + '_sorted_concat_alignments.fa', 'w') as f:
    for i in sorted_ids:
        f.write('>' + i + '\n')
        f.write(seqs[i] + '\n')