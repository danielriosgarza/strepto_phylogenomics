#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 09:30:38 2020

@author: meike
"""


'''
Testing iroki tab-seperated tree style.

'''
import matplotlib
import matplotlib.pyplot as plt

#name	branch_color	leaf_dot_color	leaf_label_color	bar1_height	bar1_gradient	bar2_height	bar2_color	new_name

def get_colors(data_list, color='coolwarm'):
    '''
    Needs a dat_list to determine how many different colours are needed. Cmap colour schemes can be provided, default are Reds.
    '''
    colors = []
    cmap = plt.get_cmap(color, len(data_list))
    
    for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))
    return colors


Tree = "(flori_1:1,(lacto_1:1,lacto_2:1,(strepto_1:1,strepto_2:1,strepto_3:0.5):0.5):0.5);"

with open ('/home/meike/files/test_tree.tre', 'w') as f:
    f.write(Tree)

#branch color = species and dot color = genus, new_name = __IROKI_BLANK__ to not be drawn

ids = ['strepto_1', 'strepto_2', 'strepto_3','lacto_1', 'lacto_2', 'flori_1']
species = ['s1', 's1','s2', 'l1', 'l2', 'f1']

unique_species = list(set(species))
species_colors = get_colors(unique_species)

species_colormapping = {}   
for i, spec in enumerate(unique_species):
    if spec not in species_colormapping:
        species_colormapping[spec] = species_colors[i]

branch_color = []
dot_color = []
new_name = []

for i, item in enumerate(ids):
    strain = species[i]
    branch_color.append(species_colormapping[strain])
    new_name.append(' __IROKI_BLANK__')
    if 'strepto' in item:
        dot_color.append('#d6604d')
    if 'lacto' in item:
        dot_color.append('#2166ac')
    if 'flori' in item:
        dot_color.append('#d1e5f0')

iroki_options = ['name', 'branch_color', 'leaf_dot_color', 'new_name']        

with open('/home/meike/files/test_tree_style.txt', 'w') as f:
    f.write('\t'.join(iroki_options) + '\n')
    for i in range(len(ids)):
        f.write(ids[i] + '\t' + branch_color[i] + '\t' + dot_color[i] + '\t __IROKI_BLANK__\n')
        
#%%

iroki_options = ['name', 'branch_color', 'leaf_label_color', 'new_name']        

with open('/home/meike/files/test_tree_style2.txt', 'w') as f:
    f.write('\t'.join(iroki_options) + '\n')
    for i in range(len(ids)):
        f.write(ids[i] + '\t' + branch_color[i] + '\t' + dot_color[i] + '\t' + species[i] + '\n')        
        