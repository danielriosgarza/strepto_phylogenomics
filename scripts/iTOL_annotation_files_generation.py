#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 08:12:20 2020

@author: meike
"""
'''
Generation of different annotation files that can be used in iTOL to visualize the phylogenetic tree.
Included: 
     >Color Range: Species colors (per genus same colour palette)
     >Outer ring marking the genus + genus labels
     >Labels: changes node names (ids) to species names
     >Bars marking genome sizes
     >Binary table showing the isolation country
     >Binary table showing the isolation source (only top 5, rest is category others)
     >Binary table of hosts
'''
from ete3 import Tree
import os
from pathlib import Path
import seaborn as sns
import random

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
ids2icountry = {}
ids2isolationsource = {}
ids2host = {}

files = [files_dir+'/03032020_streptococcus_database_final.tsv', files_dir+'/06012020_lactococcus_database.tsv', files_dir+'/06012020_floricoccus_database.tsv']


for file in files:
    with open(file) as f:
        headers = f.readline().strip().split('\t')
        species_ind = headers.index('species')
        gs_ind = headers.index('genome_length')
        icountry_ind = headers.index('isolation_country')
        isolation_ind = headers.index('isolation_source')
        host_ind = headers.index('host_name')
        for line in f:
            a = line.strip().split('\t')
            species = a[species_ind]
            ids2species[a[0]] = species
            ids2icountry[a[0]] = a[icountry_ind]
            ids2isolationsource[a[0]] = a[isolation_ind]
            ids2host[a[0]] = a[host_ind]
            if species not in species2ids:
                species2ids[species] = [a[0]]
            else:
                species2ids[species] += [a[0]]

for k, v in ids2species.items():
    id_n = int(k.split('_')[1])
    if id_n >= 11962:
        temp = v.split('-')
        species_name = 'RB-Streptococcus ' + temp[3]
        ids2species[k] = species_name

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
        
#to randomly asign sequential colors
random.shuffle(streptos)

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
# with open(annotation_files_dir +'/tree_colors_per_species.txt', 'w') as f:
#     f.write('TREE_COLORS\nSEPARATOR COMMA\nDATA\n')
#     #colored ranges inclusive labels
#     for l in leaves:
#         spe = ids2species[l.name]
#         color = leaf_colours[spe]
#         if 'strepto' in l.name:
#             f.write(l.name + ',range,' + color + ',Streptococcus\n')
#         if 'lacto' in l.name:
#             f.write(l.name + ',range,' + color + ',Lactococcus\n')
#         if 'flori' in l.name:
#             f.write(l.name + ',range,' + color + ',Floriococcus\n')
            

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
#CDS	rRNA	repeat_region	tRNA
ids2gs = {}
ids2CDS = {}
ids2rRNA = {}
ids2tRNA = {}
ids2repeat_region = {}

with open (files_dir + '/23032020_prokka_genome_data.tsv') as f:
    for line in f:
        a = line.strip().split()
        ids2gs[a[0]] = a[1]
        ids2CDS[a[0]] = a[2]
        if len(a) == 5:
            ids2rRNA[a[0]] = a[3]
            ids2repeat_region[a[0]] = a[4]
        else:
            ids2tRNA[a[0]] = a[3]

with open(annotation_files_dir + '/genomedata_RNAS_RR.txt' , 'w') as f:
    f.write('DATASET_MULTIBAR\nSEPARATOR COMMA\nDATASET_LABEL,Genome data RNAs and RR\nCOLOR,#CCE2DD\nFIELD_COLORS,#d8b365,#c7eae5,#5ab4ac,\nFIELD_LABELS,rRNA,Repeat Region,tRNA\nDATASET_SCALE,10-10-#8e9594-1-0-8,20-20-#8e9594-1-0-8,30-30-#8e9594-1-0-8,40-40-#8e9594-1-0-8,50-50-#8e9594-1-0-8,60-60-#8e9594-1-0-8,70-70-#8e9594-1-0-8,80-80-#8e9594-1-0-8,90-90-#8e9594-1-0-8,100-100-#8e9594-1-0-8\nDATA\n')
    for l in leaves:
        if l.name in ids2rRNA:
            rRNA = ids2rRNA[l.name]
            repeat_region = ids2repeat_region[l.name]
            tRNA = ''
        else:
            tRNA = ids2tRNA[l.name]
            rRNA = ''
            repeat_region = ''
        if rRNA == '':
            rRNA = 0 
        if tRNA == '':
            tRNA = 0
        if repeat_region == '':
            repeat_region = 0
        f.write(l.name + ',' + str(rRNA) + ',' + str(repeat_region) + ',' + str(tRNA) + '\n')

#Add simple bars that indicate the genome sizes at the outside of the circle
#min 927400, max 2808579
with open(annotation_files_dir + '/genomesize.txt' , 'w') as f:
    f.write('DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,Genomesizes\nCOLOR,#CCE2DD\nDATASET_SCALE,500000-500000-#8e9594-1-0-8,1000000-1000000-#8e9594-1-0-8,1500000-1500000-#8e9594-1-0-8,2000000-2000000-#8e9594-1-0-8,2500000-2500000-#8e9594-1-0-8,3000000-3000000-#8e9594-1-0-8\nDATA\n')
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

col = sns.color_palette('colorblind', n_colors = 6)
colors_shapes = col.as_hex() 
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
    if 'oral' in s or 'Oral' in s or 'subgingival' in s:
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
    elif 'gut' in s or 'Gut' in s or 'ntestine' in s or 'leostomy' in s:
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
    if len(topN) < 9:
        topN.append(count)
        indexes.append(i)
    #replace the minimum with the new count if its bigger
    else:
        min_count = min(topN)
        ind_min = topN.index(min_count)
        if count > min_count:
            topN[ind_min] = count
            indexes[ind_min] = i

#Add the top 9 sources to a list
top9 = []
for i in indexes:
    top9.append(keys[i])

source2shape = {}
for i in range(len(top9)):
    line = [str(0)]*10
    line[i] = str(1)
    source2shape[top9[i]] = ','.join(line)
source2shape['Other'] = ','.join([str(0)]*9) + ',1'

#Isolation source annotation file, only top 12 the rest is category 'others'
#9 symbols --> 1 shapes, 9 colors
slabels = list(source2shape.keys())  


col = sns.color_palette('colorblind', n_colors = len(slabels))
colors = col.as_hex() 

with open(annotation_files_dir + '/isolation_sources.txt' , 'w') as f:
    f.write('DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,Isolation Sources\nCOLOR,#01665e\nFIELD_SHAPES,2,2,2,2,2,2,2,2,2,2\nFIELD_COLORS,')
    for i, c in enumerate(colors):
        if i == colors.index(colors[-1],-1):
            f.write(c + '\nFIELD_LABELS,')
        else:
            f.write(c + ',')
    for item in slabels:
        if item == slabels[-1]:
            f.write(item + '\nDATA\n')
        else:
            f.write(item + ',')
    for l in leaves:
        original_source = ids2isolationsource[l.name]
        source = adaptedsources[original_source]
        if source in top9:
            f.write(l.name + ',' + source2shape[source] + '\n')
        elif source == '':
            pass
        else:
            f.write(l.name + ',' + source2shape['Other'] + '\n')
 
           
#Binary annotation file displaying the host names
hosts_counter = {}
for l in leaves:
    host = ids2host[l.name]
    if ',' in host:
        host = host.split(',')[1].strip()
    if host == '':
        pass
    elif host not in hosts_counter:
        hosts_counter[host] = 1
    else:
        hosts_counter[host] += 1

keys = []
values = []

for k,v in hosts_counter.items():
    keys.append(k)
    values.append(v)

#Determine topN
topN = []
indexes = []
for i, count in enumerate(values):
    #fill topN with N source counts
    if len(topN) < 9:
        topN.append(count)
        indexes.append(i)
    #replace the minimum with the new count if its bigger
    else:
        min_count = min(topN)
        ind_min = topN.index(min_count)
        if count > min_count:
            topN[ind_min] = count
            indexes[ind_min] = i

#Add the top 5 hosts to a list
top9 = []
for i in indexes:
    top9.append(keys[i])

#5 symbols: top 4 and category 'Other'
host2shape = {}
for i in range(len(top9)):
    line = [str(0)]*10
    line[i] = str(1)
    host2shape[top9[i]] = '\t'.join(line)
host2shape['Other'] = '\t'.join([str(0)]*9) + '\t1'

labels = []
for k in list(host2shape.keys()):
    if k == '': 
        k = 'Unknown'
    labels.append(k)




col = sns.color_palette('colorblind', n_colors = len(labels))
colors = col.as_hex() 

with open(annotation_files_dir + '/hosts.txt' , 'w') as f:
    f.write('DATASET_BINARY\nSEPARATOR TAB\nDATASET_LABEL\tHost names\nCOLOR\t#4d9221\nLEGEND_TITLE\tHost names\nFIELD_SHAPES\t1\t2\t3\t4\t5\t1\t2\t3\t4\t5\nFIELD_COLORS\t')
    for c in colors:
        if c == colors[-1]:
            f.write(c + '\nFIELD_LABELS\t')
        else:
            f.write(c + '\t')
    for item in labels:
        if item == labels[-1]:
            f.write(item + '\nDATA\n')
        else:
            f.write(item + '\t')
    for l in leaves:
        host = ids2host[l.name]
        if host in top9:
            f.write(l.name + '\t' + host2shape[host] + '\n')
        elif host == '':
            pass
        else:
            f.write(l.name + '\t' + host2shape['Other'] + '\n')
            
