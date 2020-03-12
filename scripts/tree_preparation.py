#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 08:12:20 2020

@author: meike
"""

from ete3 import Tree
import os
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

def get_colors(data_list, color='#d6604d'):
    '''
    Needs a dat_list to determine how many different colours are needed. Cmap colour schemes can be provided, default are Reds.
    '''
    col = sns.light_palette(color, n_colors = len(data_list))
    colors = col.as_hex()
    return colors
        
sns.set()     

path = os.getcwd()
p = Path(path)
files_dir = os.path.join(p.parents[0], 'files')

ids2species = {}

files = [files_dir+'/03032020_streptococcus_database_final.tsv', files_dir+'/06012020_lactococcus_database.tsv', files_dir+'/06012020_floricoccus_database.tsv']

for file in files:
    with open(file) as f:
        headers = f.readline().strip().split('\t')
        species_ind = headers.index('species')
        for line in f:
            a = line.strip().split('\t')
            species = a[species_ind]
            ids2species[a[0]] = species

species2ids = {}
for k, v in ids2species.items():
    species2ids[v] = k


path = os.getcwd()
p = Path(path)
files_dir = os.path.join(p.parents[0], 'files')

t = Tree('/home/meike/tests/Files/iTol/09032020_reduced_concat_alignments.fa.contree')

leaves = t.get_leaves()
ids2dups = {}
seen = {}
for leaf in leaves:
    species = ids2species[leaf.name]
    if species not in seen:
        seen[species] = leaf.name
        ids2dups[leaf.name] = []
    else:
        ids2dups[seen[species]] += [leaf.name]

with open('/home/meike/tests/Files/iTol/newcollapse.txt', 'w') as f:        
    f.write('COLLAPSE\nDATA\n')
    for k, v in ids2dups.items():
        if len(v) > 0:
            for id_ in v:
                f.write(k + '|' + id_ + '\n')

genera = {"Streptococcus" : [],
         "Lactococcus" : [],
         "Floricoccus" : []}

streptos = []
lactos = []
floris = []

for leaf in leaves:
    if 'strepto' in leaf.name:
        genera["Streptococcus"] += [leaf.name]
        streptos.append(ids2species[leaf.name])
    if 'lacto' in leaf.name:
        genera["Lactococcus"] += [leaf.name]
        lactos.append(ids2species[leaf.name])
    if 'flori' in leaf.name:
        genera["Floricoccus"] += [leaf.name]
        floris.append(ids2species[leaf.name])



strep_color = get_colors(streptos)
lacto_color = get_colors(lactos, color = '#92c5de')
flori_color = get_colors(floris, color = '#fee090')

leaf_colours = {}
for i, spec in enumerate(streptos):
    leaf_colours[spec] = strep_color[i]   

for i, spec in enumerate(lactos):
    leaf_colours[spec] = lacto_color[i]  
     
for i, spec in enumerate(floris):
    leaf_colours[spec] = flori_color[i] 


with open('/home/meike/tests/Files/iTol/tree_colors_per_species.txt', 'w') as f:
    f.write('TREE_COLORS\nSEPARATOR COMMA\nDATA\n')
    #colored ranges inclusive labels
    for l in leaves:
        spe = ids2species[l.name]
        color = leaf_colours[spe]
        if 'strepto' in l.name:
            f.write(l.name + ',range,' + color + ',Streptococcus\n')
        if 'lacto' in l.name:
            f.write(l.name + ',range,' + color + ',Lactococcus\n')
        if 'flori' in l.name:
            f.write(l.name + ',range,' + color + ',Floriococcus\n')
            



#maually check which ids are located in the middle of all species        
with open('/home/meike/tests/Files/iTol/dataset_text.txt', 'w') as f:
    #headers of the file
    f.write('DATASET_TEXT\nSEPARATOR COMMA\n')
    #set dataset labels
    f.write('DATASET_LABEL,Streptoccocus\nCOLOR,#b2182b\nALIGN_TO_TREE,1\n')
    f.write('DATASET_LABEL,Lactococcus\nCOLOR,#2166ac\nALIGN_TO_TREE,1\n')
    f.write('DATASET_LABEL,Floricoccus\nCOLOR,#525252\nALIGN_TO_TREE,1\n')
    f.write('DATA\n')
    f.write('streptococcus_11897,Streptococcus,-1,#b2182b,bold-italic,3,0\n')
    f.write('lactococcus_00017,Lactococcus,-1,#2166ac,bold-italic,3,0\n')
    f.write('floricoccus_00001,Floricoccus,-1,#525252,bold-italic,3,0\n')
    

#look at tree and search for ids from each genus that are furthest apart
all_species = []
for leaf in leaves:
    spec = ids2species[leaf]

                   
with open('/home/meike/tests/Files/iTol/tree_colors.txt', 'w') as f:
    f.write('TREE_COLORS\nSEPARATOR COMMA\nDATA\n')
    #colored ranges inclusive labels
    f.write('streptococcus_11918|streptococcus_01671,range,#f9665e,Streptococcus\n')
    f.write('lactococcus_00137|lactococcus_00009,range,#799fcb,Lactococcus\n')
    f.write('floricoccus_00001,range,#eef1e6,Floricoccus\n')

#Change labels of nodes    
with open('/home/meike/tests/Files/iTol/labels.txt', 'w') as f:
    f.write('LABELS\nSEPARATOR COMMA\nDATA\n')
    for k,v in genera.items():
        for item in v:
            species = ids2species[item]
            f.write(item + ',' + species + '\n')
            
with open('/home/meike/tests/Files/iTol/dataset_style.txt', 'w') as f:
    f.write('DATASET_STYLE\nSEPARATOR COMMA\n')
    f.write('DATASET_LABEL,Streptoccocus\nCOLOR,#f9665e\n')
    f.write('DATASET_LABEL,Lactococcus\nCOLOR,#799fcb\n')
    f.write('DATASET_LABEL,Floricoccus\nCOLOR,#eef1e6\n')
    f.write('LEGEND_TITLE,Genus\nLEGEND_SHAPES,2,2,2\nLEGEND_COLORS,#f9665e,#799fcb,#eef1e6\nLEGEND_LABELS,Streptoccocus,Lactococcus,Floricoccus\nDATA\n')
    for k,v in genera.items():
        for item in v:
            if k == "Streptococcus":
                f.write(item + ',label,node,#000000,1,normal,#f9665e\n')
            if k == "Lactococcus":
                f.write(item + ',label,node,#000000,1,normal,#799fcb\n')
            if k == "Floricoccus":
                f.write(item + ',label,node,#000000,1,normal,#eef1e6\n')
                



