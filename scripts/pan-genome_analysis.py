#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 13:35:48 2020

@author: meike
"""

'''
Relation of Pan-genome and Phylogenomics.
Step1:
    >Grouping the Strepto_ids (8 according to distances in tree)
    >Done manually: open tree in iTOL and copy the ids belonging to a certain ancestor node to clipboard. Save all in a file and make table:
        all ids in rows and group 1-8 as columns. 1 = part of this group, 0 = absent in this group.
        
Step2:
    >Fishers exact test
    
Step3:
    >Heatmap
    
'''

import os
from pathlib import Path
from datetime import date
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import scipy.stats as sts
import seaborn as sns

def subsample_fisher_exact(gene_strain, gene_out_strain, n_samples):
    samp_size = len(gene_strain)
    pvalues_greater = np.zeros(n_samples)
    pvalues_less = np.zeros(n_samples)
    
    for i in range(n_samples):
        rchoice= np.random.choice(np.arange(len(gene_out_strain)), size=samp_size)
        samp = gene_out_strain[rchoice]
        pvalues_greater[i] = min(1,sts.fisher_exact([[sum(gene_strain), len(gene_strain)-sum(gene_strain)], [sum(samp), len(samp)-sum(samp)]], alternative='greater')[1]*2)
        pvalues_less[i] = min(sts.fisher_exact([[sum(gene_strain), len(gene_strain)-sum(gene_strain)], [sum(samp), len(samp)-sum(samp)]], alternative = 'less')[1]*2,1)
    
    
    return sts.hmean(pvalues_greater),sts.hmean(pvalues_less)


path = os.getcwd()
p = Path(path)

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)

#%% runcell 1
#Step 1: defining Streptococcus subgroups


#make dict of ids mapping to its group, and a dict of groups mapping to binary 'code' that is used later for the table
groups = {}
bin_table = {}
num_group = 0

with open(os.path.join(p.parents[0], 'files', 'binary_table', 'groups_tree_streptococcus.txt')) as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            group = num_group
            code = ['0']*8
            code[num_group] = str(num_group)
            bin_table[group] = '\t'.join(code)
            num_group += 1
        else:
            a = line.split(' ')
            id_ = '_'.join(a)
            groups[id_] = group

#write ids into table (rows) and mark the according group with a 1 (columns)

ids = sorted(list(groups.keys()))

with open(os.path.join(p.parents[0], 'files', 'binary_table', today + 'strepto_groups.tsv'), 'w') as f:
    f.write('genome_id\tgroup1\tgroup2\tgroup3\tgroup4\tgroup5\tgroup6\tgroup7\tgroup8\n')
    for i in ids:
        group = groups[i]
        f.write(i + '\t' + bin_table[group] + '\n')
        
        
#%% runcell 2

#Step 2: 
#Exploratory testing: test all genes that are at least as much present in a number of strains as the smallest group (e.g. smallest group = 10, genes should be at least present in 10 strains)
# Use fishers exact test: compare gene present in group to not present in group, present in rest, not present in rest (rest: random choices of all other groups (size: as big as the group that is tested))

bin_file = os.path.join(p.parents[0], 'files', 'binary_table', '14042020_binary_table_sorted.tsv')

#determine smallest group 
counts = {}
for k, v in groups.items():
    if v not in counts:
        counts[v] = 1
    else:
        counts[v] += 1
#print(counts) #smallest group is group 2 with 10 strains

groups_array = np.zeros(len(ids))
for i, id_ in enumerate(ids):
    groups_array[i] = groups[id_]

#look at gene presence (binary part of sorted table)
with open(bin_file) as f:
    lines = [line for line in f]
ncols = len(lines[0].split('\t'))

data = np.loadtxt(bin_file, delimiter = '\t', skiprows = 1, usecols = range(5, ncols))

#subset of genes that is at least present in 10 strains (equal to smallest group) used as cutoff to exclude artefical significance
genes_and_strains = data[np.where(np.count_nonzero(data, axis = 1) >= 10)]

#test_genes = genes_and_strains[np.random.randint(genes_and_strains.shape[0], size = 5), :]

genes = [g for g in range(len(genes_and_strains))]

#make p-value Dict to save all pvalues per group per gene 
pv_dict={}

for i,gr in enumerate(groups_array):
    if gr not in pv_dict:   
        pv_dict[gr]={}
    for idx_ge, g in enumerate(genes):
        gene_pab = genes_and_strains[idx_ge]
        gene_strain = gene_pab[groups_array == gr]
        gene_out_strain = gene_pab[groups_array != gr]
        
        #fisher exact: group vs random choice of rest of groups, these random choices are done 10 times and the hmean of the pvalue is used and added to pv_dict
        pv_dict[gr][g] = subsample_fisher_exact(gene_strain, gene_out_strain, 10) 
 
        
#write two files with pvalues of genes: one with pvalues indicating greater presence and one with less presence in one group compared to others
#gene order is the same as in sorted binary table

group_order = sorted(list(pv_dict.keys()))

with open(os.path.join(p.parents[0], 'files', 'binary_table', today + '_pvalues_genes_greater.tsv'), 'w') as f:
    f.write('gene_number\tgroup1\tgroup2\tgroup3\tgroup4\tgroup5\tgroup6\tgroup7\tgroup8\n')
    for ge in genes:
        line = 'gene_' + str(ge) + '\t'
        for gr in group_order:
            
                line += '%.5f' % pv_dict[gr][ge][0] + '\t'
        f.write(line + '\n')
        
with open(os.path.join(p.parents[0], 'files', 'binary_table', today + '_pvalues_genes_less.tsv'), 'w') as f:
    f.write('gene_number\tgroup1\tgroup2\tgroup3\tgroup4\tgroup5\tgroup6\tgroup7\tgroup8\n')
    for ge in genes:
        line = 'gene_' + str(ge) + '\t'
        for gr in group_order:
            line += '%.5f' % pv_dict[gr][ge][1] + '\t'
        f.write(line + '\n')


#%% runcell 3

#step 3: Make heat map

#make heatmap from -log 10 p-values to indicate significant differences between groups (greater = gene more present in this group compared to rest, less = gene less present in groupo compared to rest)

pvs_greater = {}
pvs_less = {}

for gr in pv_dict:
    pvs_greater[gr] = {}
    pvs_less[gr] = {}
    for ge in pv_dict[gr]:
        pvs_greater[gr][ge] = -np.log10(pv_dict[gr][ge][0])
        pvs_less[gr][ge] = -np.log10(pv_dict[gr][ge][1])

data_vis_greater = np.zeros((len(pvs_greater), len(genes)))
data_vis_less = np.zeros((len(pvs_less), len(genes)))

#to keep the groups in the right order use a sorted list of the dict keys
for i_gr, gr in enumerate(group_order):
    for i_ge, ge in enumerate(genes):
        data_vis_greater[i_gr][i_ge] = pvs_greater[gr][ge]
        data_vis_less[i_gr][i_ge] = pvs_less[gr][ge]

sum_rows = np.sum(data_vis_greater,axis=0)
sum_cols = np.sum(data_vis_greater,axis=1)
sort_rows = np.argsort(sum_rows)
sort_rows=sort_rows[::-1]
sort_cols = np.argsort(sum_cols)
sort_cols=sort_cols[::-1]

l_sum_rows = np.sum(data_vis_less,axis=0)
l_sum_cols = np.sum(data_vis_less,axis=1)
l_sort_rows = np.argsort(l_sum_rows)
l_sort_rows = l_sort_rows[::-1]
l_sort_cols = np.argsort(l_sum_cols)
l_sort_cols = l_sort_cols[::-1]

#Heatmap of genes more present in certain group
fig1, ax1 = plt.subplots(1, figsize= (10,10))

ax1 = sns.heatmap(data_vis_greater[sort_cols].T[sort_rows], cmap = 'coolwarm', linewidths=0.1, linecolor='w')

xlab = list(np.array([0,1,2,3,4,5,6,7])[sort_cols])
ax1.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5], xlab)
ax1.set_xticklabels(xlab)

ax1.set_title('Genes with greater presence in certain group\n(-log10(p-value))')


#heatmap of genes less present compared to other groups
fig2, ax2 = plt.subplots(1, figsize= (10,10))
ax2 = sns.heatmap(data_vis_less[l_sort_cols].T[l_sort_rows], cmap = 'coolwarm', linewidths=0.1, linecolor='w')

l_xlab = list(np.array([0,1,2,3,4,5,6,7])[l_sort_cols])
ax2.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5], l_xlab)
ax2.set_xticklabels(l_xlab)

ax2.set_title('Genes with less presence in certain group\n(-log10(p-value))')

fig1.tight_layout()
fig1.savefig(os.path.join(p.parents[0], 'figures', today + '_heatmap_greater_genes.png'), dpi=300, bbox_inches='tight')

fig2.tight_layout()
fig2.savefig(os.path.join(p.parents[0], 'figures', today + '_heatmap_less_genes.png'), dpi=300, bbox_inches='tight')
