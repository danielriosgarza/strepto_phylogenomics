#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 09:25:36 2020

@author: meike
"""

'''
Look per group which genes have p_value >= 0.00005. Get binary lines. Make binary of genes to groups: 8 groups x 2 directions (greater/less) = 16 total
'''

import os
from pathlib import Path
from datetime import date
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def binplot(Array, Title):
    '''
    makes binary plots. Gene present = light blue(1), absent = black(0)
    '''
    fig, ax = plt.subplots()

    # define the colors
    cmap = mpl.colors.ListedColormap(['#404040', '#d1e5f0'])


    fig.suptitle(Title)
    # plot it
    ax.imshow(array, interpolation='none', cmap=cmap)



path = os.getcwd()
p = Path(path)

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)


#%% runcell 1

#Looking for genes that are specifically for a certain group
#Make binary plots and save specific genes in lists to write them into files (runcell 2) and into summary tables (runcell 3)

#Needed files
bin_file = os.path.join(p.parents[0], 'files', 'binary_table', '14042020_binary_table_sorted.tsv')
    
bin_greater = os.path.join(p.parents[0], 'files', 'binary_table', '24042020_pvalues_genes_greater.tsv')

bin_less = os.path.join(p.parents[0], 'files', 'binary_table', '24042020_pvalues_genes_less.tsv')

group_file = os.path.join(p.parents[0], 'files', 'binary_table','04052020strepto_groups.tsv')


#Dict mapping each id to group it belongs to                          
ids2group = {}  
with open (group_file) as f:
        f.readline()
        for line in f:
            a = line.strip().split('\t')
            for item in a[1::]:
                if item != '0':
                    ids2group[a[0]] = item

#and other way around
group2ids = {}
for k, v in ids2group.items():
    if v not in group2ids:
        group2ids[v] = [k]
    else:
        group2ids[v].append(k)
    
#Dict with all groups mapping to gene numbers (!) with greater presence --> look for gene names in sorted binary table for gene lists
group2grgenes = {}
with open(bin_greater) as f:
    f.readline()
    for line in f:
        a = line.strip().split('\t')
        for ind, value in enumerate(a[1::]):
            groupnum = ind + 1
            genenum = int(a[0].split('_')[1])
            
            #used a cutoff of 5x10-5 to get list of genes that is probably specific for a group
            if float(value) <= 0.00005:
                if groupnum not in group2grgenes:
                    group2grgenes[groupnum] = [genenum]
                else:
                    group2grgenes[groupnum] = group2grgenes[groupnum] + [genenum]
  


#get binary lines to plot and save gene lists

with open(bin_file) as f:
    ids = f.readline().strip().split('\t')[5::]
    lines = [line.strip().split('\t') for line in f]

group2unigenes = {}    
group2binlines = {}
for gr in group2grgenes:
    group2binlines[gr] = []
    group2unigenes[gr] = []
    
    #Go through all gene indexes to get binary lines that show, which id has this gene
    for ind in group2grgenes[gr]:
        
        #make binary for presence of gene in groups 1-8
        groups = [0]*8 
        
        #Look through the binary of all strains if they contain the gene (lines[ind] = complete line at position of the gene)
        for i, value in enumerate(lines[ind][5::]):
            id_ = ids[i]
            gr_num = ids2group[id_]
            if value == '1':
                
                #if gene is present in strain of a certain group mark it in binary with a 1 (needs correct index: group 1 = index 0)
                groups[int(gr_num) - 1] = 1
                
        #look if gene is only present in the group you look at
        if groups.count(1) == 1 and ind not in group2unigenes[gr]:
            group2unigenes[gr].append(ind)
                
        #add binary line to Dict
        group2binlines[gr].append(groups)
        
        #add group to Dicts if it has no specific genes
        for n in '12345678':
            if int(n) not in group2unigenes:
                group2unigenes[int(n)] = []
                group2binlines[int(n)] = [0]*8
                
    #make array out of the binary lines to plot them
    array = np.array(group2binlines[gr], dtype= 'float')
    binplot(array, str(gr) + ' greater')


#same things only for less present genes

lgroup2grgenes = {}
with open(bin_less) as f:
    f.readline()
    for line in f:
        a = line.strip().split('\t')
        for ind, value in enumerate(a[1::]):
            groupnum = ind + 1
            genenum = int(a[0].split('_')[1])
            if float(value) <= 0.00005:
                if groupnum not in lgroup2grgenes:
                    lgroup2grgenes[groupnum] = [genenum]
                else:
                    lgroup2grgenes[groupnum] = lgroup2grgenes[groupnum] + [genenum]
  

lgroup2unigenes = {}    
lgroup2binlines = {}
for gr in lgroup2grgenes:
    lgroup2binlines[gr] = []
    lgroup2unigenes[gr] = []
    for ind in lgroup2grgenes[gr]:
        #binary for presence of gene in groups 1-8, binary other way around: should not be present in group of interest (less presence of gene)
        groups = [1]*8 
        for i, value in enumerate(lines[ind][5::]):
            id_ = ids[i]
            gr_num = ids2group[id_]
            if value == '0':
                
                #if gene is present in strain of a certain group mark it in binary with a 0 (needs correct index: group 1 = index 0)
                groups[int(gr_num) - 1] = 0
        if groups.count(0) == 1 and ind not in lgroup2unigenes[gr]:
                lgroup2unigenes[gr].append(ind)
        lgroup2binlines[gr].append(groups)


    array = np.array(lgroup2binlines[gr], dtype= 'float')
    binplot(array, str(gr) + ' less')

#%% runcell 2

#Save lists of unique genes per group in a file

with open(os.path.join(p.parents[0], 'files', today + '_unique_genes_per_group.tsv'), 'w') as f:
    f.write('Group_number\tgreater_presence\tless_presence\n')
   
    #go through list of greater genes per group (only numbers yet) 
    for gr in sorted(list(group2unigenes.keys())):
        f.write('group_' + str(gr) + '\t')
        
        #get gene names and if not available, a protein description
        gene_descrip = []
        prot_d = []
        
        #if list of unique genes empty = no specific genes: add N/A
        if len(group2unigenes[gr]) == 0:
            f.write('N/A\t')
        else:
            for ind in group2unigenes[gr]:
                gene = lines[ind][2]
                prot_info = lines[ind][4]
                
                #only add gene names once and not N/A
                if gene != 'N/A' and gene not in gene_descrip:
                    gene_descrip.append(gene)
                if prot_info not in prot_d:
                    prot_d.append(prot_info)
                    
            #Check if gene name were found
            if len(gene_descrip) >= 1:
                f.write(','.join(gene_descrip) + '\t')
            else:
                f.write('protein_description: ' + ','.join(prot_d) + '\t')
                
                
        #same for less present
        lgene_descrip = []
        lprot_d = []
        if len(lgroup2unigenes[gr]) == 0:
            f.write('N/A\n')
        else:
            for ind in lgroup2unigenes[gr]:
                lgene = lines[ind][2]
                lprot_info = lines[ind][4]
                if lgene != 'N/A' and lgene not in lgene_descrip:
                    lgene_descrip.append(lgene)
                if lprot_info not in lprot_d:
                    lprot_d.append(lprot_info)
            if len(lgene_descrip) >= 1:
                f.write(','.join(lgene_descrip) + '\n')
            else:
                f.write('protein_description: ' + ','.join(lprot_d) + '\n')
                
#%% runcell 3

#Summarizing table for all groups

#Make Dict with all ids mapping to the species of the genome
ids2species = {}

files_dir = os.path.join(p.parents[0], 'files')
files = [files_dir+'/03032020_streptococcus_database_final.tsv', files_dir+'/06012020_lactococcus_database.tsv', files_dir+'/06012020_floricoccus_database.tsv']


for file in files:
    with open(file) as f:
        headers = f.readline().strip().split('\t')
        species_ind = headers.index('species')

        for line in f:
            a = line.strip().split('\t')
            species = a[species_ind]
            ids2species[a[0]] = species

#Change names of the ones from Nijmegen and uncultured --> unclassified            
for k, v in ids2species.items():
    id_n = int(k.split('_')[1])
    if  'uncultured' in v:
        ids2species[k] = 'Unclassified Streptococcus'
                
    if id_n >= 11962:
        temp = v.split('-')
        species_name = 'Streptococcus ' + temp[3]
        ids2species[k] = species_name

    
#get gene names per group
with open(os.path.join(p.parents[0], 'files', today + '_unique_genes_per_group.tsv')) as f:
    f.readline()
    gene_lines = [l.strip().split('\t')[1::] for l in f]

#Make list of all species to set it later in alphabeticl order
all_specs = []
s2gr = {}
for i in ids:
    s = ids2species[i]
    gr = ids2group[i]
    if s not in s2gr:
        s2gr[s] = gr
    if s not in all_specs:
        all_specs.append(s)
all_specs = sorted(all_specs)

#write summary file
with open(os.path.join(p.parents[0], 'files', today + '_summary_groups.tsv'), 'w') as f:
    f.write('Species\tGroup\tGenes with greater presence\tGenes with less presence\n')
    for sp in all_specs:
        gr_n = s2gr[sp]
        f.write(sp + '\tGroup_' + gr_n + '\t' +  gene_lines[int(gr_n)-1][0] + '\t' + gene_lines[int(gr_n)-1][1] + '\n')
        
##Other order of summary table:
# with open(os.path.join(p.parents[0], 'files', today + '_summary_groups.tsv'), 'w') as f:
#     f.write('Group\tSpecies\tGenes with greater presence\tGenes with less presence\n')
#     for n in '12345678':
        
#         #First get species belonging to the group in alphabetical order
#         specs = []
#         for id_ in group2ids[n]:
#             spec = ids2species[id_]
#             if spec not in specs:
#                 specs.append(spec)
#         specs = sorted(specs)
        
#         #write line: group number, species, genes greater, genes less
#         f.write('Group ' + n + '\t' + ','.join(specs) + '\t' +  gene_lines[int(n)-1][0] + '\t' + gene_lines[int(n)-1][1] + '\n')
