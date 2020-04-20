#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 09:58:59 2020

@author: daniel
"""

from pylab import *
import scipy.stats as sts
import numpy as np
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
    


genes_and_strains = np.array([[0., 0., 0., 1., 1., 1.], [0., 0., 1., 1., 1., 0.], [1., 0., 0., 1., 1., 0.], [1., 1., 1., 0., 1., 1.], [0., 1., 0., 1., 1., 1.]])

groups = np.array([0,0,0,1,1,1])

genes = np.arange(5)


pv_dict={}

for i,gr in enumerate(groups):
    pv_dict[gr]={}
    for idx_ge, g in enumerate(genes):
        gene_pab = genes_and_strains[idx_ge]
        gene_strain = gene_pab[groups==gr]
        gene_out_strain = gene_pab[groups!=gr]
        
        
        
        pv_dict[gr][g] = subsample_fisher_exact(gene_strain, gene_out_strain,10)
        

single_v_pv={}

for gr in pv_dict:
    single_v_pv[gr]={}
    for ge in pv_dict[gr]:
        single_v_pv[gr][ge] = -np.log10(pv_dict[gr][ge][0])


data = np.zeros((len(single_v_pv), len(genes)))

for i_gr, gr in enumerate(pv_dict):
    for i_ge, ge in enumerate(genes):
        data[i_gr][i_ge] = single_v_pv[gr][ge]

sum_rows = np.sum(data,axis=0)
sum_cols = np.sum(data,axis=1)
sort_rows = np.argsort(sum_rows)
sort_rows=sort_rows[::-1]
sort_cols = np.argsort(sum_cols)
sort_cols=sort_cols[::-1]

sns.heatmap(data[sort_cols].T[sort_rows], cmap=cm.coolwarm, linewidths=0.1, linecolor='w')

xlab = np.array([0,1])[sort_cols]
xticks([0.5,1.5], xlab)