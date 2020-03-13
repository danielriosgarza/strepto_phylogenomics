#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 08:12:20 2020

@author: meike
"""

from ete3 import Tree
import os
from pathlib import Path
import seaborn as sns

def get_colors(data_list, color='#d6604d'):
    '''
    Needs a data_list to determine how many different colours are needed. Hexcode color can be provided, default is a light red.
    '''
    col = sns.light_palette(color, n_colors = len(data_list))
    colors = col.as_hex()
    return colors
        
sns.set()     

path = os.getcwd()
p = Path(path)
files_dir = os.path.join(p.parents[0], 'files')


#Make Dict with all ids mapping to the species of the genome, and a Dict where species are mapping to ids
ids2species = {}
species2ids = {}
files = [files_dir+'/03032020_streptococcus_database_final.tsv', files_dir+'/06012020_lactococcus_database.tsv', files_dir+'/06012020_floricoccus_database.tsv']

for file in files:
    with open(file) as f:
        headers = f.readline().strip().split('\t')
        species_ind = headers.index('species')
        for line in f:
            a = line.strip().split('\t')
            species = a[species_ind]
            ids2species[a[0]] = species
            if species not in species2ids:
                species2ids[species] = [a[0]]
            else:
                species2ids[species] += [a[0]]

#open the tree file and get information over included ids (leaves)
t = Tree('/home/meike/tests/Files/iTol/09032020_reduced_concat_alignments.fa.contree')

leaves = t.get_leaves()


#Get colors of from same palette for each genus
streptos = []
lactos = []
floris = []

strep_color = get_colors(streptos)
lacto_color = get_colors(lactos, color = '#92c5de')
flori_color = get_colors(floris, color = '#969696')

leaf_colours = {}
for i, spec in enumerate(streptos):
    leaf_colours[spec] = strep_color[i]   

for i, spec in enumerate(lactos):
    leaf_colours[spec] = lacto_color[i]  
     
for i, spec in enumerate(floris):
    leaf_colours[spec] = flori_color[i] 

#write annotation file for iTOL: colors background of the leaves according to its species
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
            

#Label annotation file iTOL: manually check where the label should be located and choose id (leaf) that should be labled        
with open('/home/meike/tests/Files/iTol/dataset_text.txt', 'w') as f:
    #headers of the file
    f.write('DATASET_TEXT\nSEPARATOR COMMA\n')
    #set dataset labels
    f.write('DATASET_LABEL,Genus labels\nCOLOR,#b2182b\nALIGN_TO_TREE,1\n')
    f.write('DATA\n')
    f.write('streptococcus_11897,Streptococcus,-1,#c23b22,bold-italic,3,0\n')
    f.write('lactococcus_00017,Lactococcus,-1,#2166ac,bold-italic,3,0\n')
    f.write('floricoccus_00001,Floricoccus,-1,#969696,bold-italic,3,0\n')
    


#label annotation file: changes the labels at leaf_nodes   
with open('/home/meike/tests/Files/iTol/labels.txt', 'w') as f:
    f.write('LABELS\nSEPARATOR COMMA\nDATA\n')
    for l in leaves:
        species = ids2species[l.name]
        f.write(l.name + ',' + species + '\n')


#Colorstrip (outer ring) annotation file: strip at the outside is colored according to the genus a node belongs to                
with open('/home/meike/tests/Files/iTol/dataset_colorstrip.txt', 'w') as f:
    f.write('DATASET_COLORSTRIP\nSEPARATOR COMMA\n')
    f.write('DATASET_LABEL,Species Colors\nCOLOR,#f9665e\n')
    f.write('COLOR_BRANCHES,1\nDATA\n')
    for l in leaves:
        if 'strepto' in l.name:
            f.write(l.name + ',#c23b22,Streptococcus\n')
        if 'lacto' in l.name:
            f.write(l.name + ',#2166ac,Lactococcus\n')
        if 'flori' in l.name:
            f.write(l.name + ',#969696,Floriococcus\n')

