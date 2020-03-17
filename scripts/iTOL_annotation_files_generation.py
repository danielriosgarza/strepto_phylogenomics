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
annotation_files_dir = os.path.join(p.parents[0], 'files', 'phylogenetic_tree', 'iTOL')

#Make Dict with all ids mapping to the species of the genome, and a Dict where species are mapping to ids
ids2species = {}
species2ids = {}
ids2gs = {}
ids2icountry = {}
ids2isolationsource = {}

files = [files_dir+'/03032020_streptococcus_database_final.tsv', files_dir+'/06012020_lactococcus_database.tsv', files_dir+'/06012020_floricoccus_database.tsv']


for file in files:
    with open(file) as f:
        headers = f.readline().strip().split('\t')
        species_ind = headers.index('species')
        gs_ind = headers.index('genome_length')
        icountry_ind = headers.index('isolation_country')
        isolation_ind = headers.index('isolation_source')
        for line in f:
            a = line.strip().split('\t')
            species = a[species_ind]
            ids2species[a[0]] = species
            ids2gs[a[0]] = a[gs_ind]
            ids2icountry[a[0]] = a[icountry_ind]
            ids2isolationsource[a[0]] = a[isolation_ind]
            if species not in species2ids:
                species2ids[species] = [a[0]]
            else:
                species2ids[species] += [a[0]]

#open the tree file and get information over included ids (leaves)
t = Tree(files_dir + '/phylogenetic_tree/12032020_reduced_concat_alignments.fa.contree')

leaves = t.get_leaves()
    
#Get colors of from same palette for each genus
streptos = []
lactos = []
floris = []

for l in leaves:
    v = ids2species[l.name]
    if 'strepto' in l.name and v not in streptos:
        streptos.append(v)
    if 'lacto' in l.name and v not in lactos:
        lactos.append(v)
    if 'flori' in l.name and v not in floris:
        floris.append(v)

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
with open(annotation_files_dir +'/tree_colors_per_species.txt', 'w') as f:
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
with open(annotation_files_dir + '/dataset_text.txt', 'w') as f:
    #headers of the file
    f.write('DATASET_TEXT\nSEPARATOR COMMA\n')
    #set dataset labels
    f.write('DATASET_LABEL,Genus labels\nCOLOR,#b2182b\nALIGN_TO_TREE,1\n')
    f.write('DATA\n')
    f.write('streptococcus_11897,Streptococcus,-1,#c23b22,bold-italic,3,0\n')
    f.write('lactococcus_00183,Lactococcus,-1,#2166ac,bold-italic,3,0\n')
    f.write('floricoccus_00001,Floricoccus,-1,#969696,bold-italic,3,0\n')
    


#label annotation file: changes the labels at leaf_nodes   
with open(annotation_files_dir + '/labels.txt', 'w') as f:
    f.write('LABELS\nSEPARATOR COMMA\nDATA\n')
    for l in leaves:
        species = ids2species[l.name]
        f.write(l.name + ',' + species + '\n')


#Colorstrip (outer ring) annotation file: strip at the outside is colored according to the genus a node belongs to                
with open(annotation_files_dir + '/dataset_colorstrip.txt', 'w') as f:
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

#Add simple bars that indicate the genome sizes at the outside of the circle
#min 927400, max 2808579
with open(annotation_files_dir + '/genomesize_bars.txt' , 'w') as f:
    f.write('DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,Genomesizes\nCOLOR,#CCE2DD\nDATASET_SCALE,0,500000,1000000,1500000,2000000,2500000,3000000\nDATA\n')
    for l in leaves:
        gs = ids2gs[l.name]
        if gs == '':
            gs = 0 
        f.write(l.name + ',' + str(gs) + '\n')
        
#Binary dataset for isolation country
countries = []        
for l in leaves:
    ic = ids2icountry[l.name]
    countries.append(ic)
    
ndifferent_shapes = len(set(countries)) #25 minus empty spaces = 24

#24 different legend symbols --> 4 shapes x 6 different colors

colors_shapes = ['#b2182b', '#92c5de', '#d6604d', '#2166ac', '#f4a582', '#4393c3']
labels = list(set(countries))
labels.remove('')

country2shape = {}
for i in range(len(labels)):
    line = [str(0)]*24
    line[i] = str(1)
    country2shape[labels[i]] = ','.join(line)
country2shape[''] = ','.join([str(0)]*24)
      
with open(annotation_files_dir + '/isolation_countries.txt' , 'w') as f:
    f.write('DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,Isolation Countries\nCOLOR,#B05F3C\nFIELD_SHAPES,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4\nFIELD_COLORS,')
    for i in range(6):
        for j in range(6):
            if i == 5 and j == 5:
                f.write(colors_shapes[j] + '\n')
            else:
                f.write(colors_shapes[j] + ',')
    f.write('FIELD_LABELS,')
    for item in labels:
        if item == labels[-1]:
            f.write(item + '\n')
        else:
            f.write(item + ',')
    f.write('DATA\n')
    for l in leaves:
        c = ids2icountry[l.name]
        if c != '':
            f.write(l.name + ',' + country2shape[c] + '\n')

#Binary dataset for isolation sources. Only display top 5 and rest in category 'others'
            
sources = []
for l in leaves:
    sources.append(ids2isolationsource[l.name])

adaptedsources = {}
sources_counter = {}

for s in sources:
    if 'oral' in s or 'Oral' in s:
        adaptedsources[s] = 'Oral'
        if 'Oral' not in sources_counter:
            sources_counter['Oral'] = 1
        else:
            sources_counter['Oral'] += 1
    elif 'airy' in s or 'ilk' in s:
        adaptedsources[s] = 'Dairy'
        if 'Dairy' not in sources_counter:
            sources_counter['Dairy'] = 1
        else:
            sources_counter['Dairy'] += 1
    elif 'lood' in s:
        adaptedsources[s] = 'Blood'
        if 'Blood' not in sources_counter:
            sources_counter['Blood'] = 1
        else:
            sources_counter['Blood'] += 1
    elif 'gut' in s or 'Gut' in s or 'ntestine' in s:
        adaptedsources[s] = 'Intestine'
        if 'Intestine' not in sources_counter:
            sources_counter['Intestine'] = 1
        else:
            sources_counter['Intestine'] += 1
    elif 'feces' in s or 'fecal' in s or 'stool' in s:
        adaptedsources[s] = 'Feces'
        if 'Feces' not in sources_counter:
            sources_counter['Feces'] = 1
        else:
            sources_counter['Feces'] += 1
    else:
        adaptedsources[s] = s.capitalize()
        if s.capitalize() not in sources_counter:
            sources_counter[s.capitalize()] = 1
        else:
            sources_counter[s.capitalize()] += 1

del sources_counter['']

keys = []
values = []

for k,v in sources_counter.items():
    keys.append(k)
    values.append(v)

#Determine topN
topN = []
indexes = []
for i, count in enumerate(values):
    #fill topN with N source counts
    if len(topN) < 5:
        topN.append(count)
        indexes.append(i)
    #replace the minimum with the new count if its bigger
    else:
        min_count = min(topN)
        ind_min = topN.index(min_count)
        if count > min_count:
            topN[ind_min] = count
            indexes[ind_min] = i

#Add the top 5 sources to a list
top5 = []
for i in indexes:
    top5.append(keys[i])

source2shape = {}
for i in range(len(top5)):
    line = [str(0)]*6
    line[i] = str(1)
    source2shape[top5[i]] = ','.join(line)
source2shape['Other'] = ','.join([str(0)]*5) + ',1'
#Isolation source annotation file, only top 5 the rest is category 'others'
#6 symbols --> 6 shapes
  
labels = list(source2shape.keys())  
with open(annotation_files_dir + '/isolation_sources.txt' , 'w') as f:
    f.write('DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,Isolation Sources\nCOLOR,#01665e\nFIELD_SHAPES,1,2,3,4,5,2\nFIELD_COLORS,#3182bd,#3182bd,#3182bd,#3182bd,#3182bd,#bd0026\nFIELD_LABELS,')
    for item in labels:
        if item == labels[-1]:
            f.write(item + '\nDATA\n')
        else:
            f.write(item + ',')
    for l in leaves:
        original_source = ids2isolationsource[l.name]
        source = adaptedsources[original_source]
        if source in top5:
            f.write(l.name + ',' + source2shape[source] + '\n')
        else:
            f.write(l.name + ',' + source2shape['Other'] + '\n')
        