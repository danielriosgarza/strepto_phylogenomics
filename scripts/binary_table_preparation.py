#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 11:55:25 2020

@author: meike
"""


'''
Binary table of Pan-genome. Product information of all orthologues found within a group are stated in the first column. The information is obtained from the prokka annotation (<id>.tsv) Other columns are genome ids displaying if certain protein is present (1) or absent (0). 
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


taxon_list = []
with open(os.path.join(p.parents[0], 'files', 'taxon_list')) as f:
    for line in f:
        taxon_list.append(line.strip())

        
#make binary table: all orthologs (product info) in a list with all genomes with 0 and 1 following
   
with open (os.path.join(p.parents[0], 'files', 'binary_table', 'all.ort.group')) as f:
    with open(os.path.join(p.parents[0], 'files', 'binary_table', today + '_binary_table_prep.tsv'), 'w') as f2:
        
        #header: write orthologues followed by all genome ids for the binary table
        f2.write('Orthologue\t')
        
        for db_id in taxon_list:
            if db_id == taxon_list[-1]:
                f2.write(db_id + '\n')
            else:
                f2.write(db_id + '\t')
                
        #go through mcl file and look through each orthogroup (=line)
        for line in f:
            a = line.strip().split('\t')
            ids = []
            orthos = []
            
            #go through all ids and get product info of the orthologue from the prokka annotation file    
            for pair in a:
                id_ = pair.split('|')[0]
                ids.append(id_)
                ortho = pair.split('|')[1]
                orthos.append(ortho)
            
            #reduce the amount of unique and cloud genes (present in less than 0.05%)
            if len(set(ids)) > 60:
                #make dict with all ids mapping to all orthologues/paralogues to get protein information
                ids2orthos = {}
                for i, num in enumerate(ids):
                    if num not in ids2orthos:
                        ids2orthos[num] = [orthos[i]]
                    else:
                        ids2orthos[num] += [orthos[i]]
                #get product information from the prokka .tsv file
                products = []        
                for k, v in ids2orthos.items():        
                    with open('/home/meiker/git/data/prokka_annotation/' + k + '/' + k + '.tsv') as f3:
                        f3.readline()
                        for line in f3:
                            a = line.strip().split('\t')
                            if a[0] in v and a[-1] not in products:
                                products.append(a[-1])
                
                #list first all information of protein function
                f2.write(','.join(products) + '\t')
                
                #mark found ids with an 1 in the table, rest is 0 
                group = [str(0)]*len(taxon_list)
                for id_ in set(ids):
                    group[taxon_list.index(id_)] = str(1)
                f2.write('\t'.join(group) + '\n')
                
            
                
            
        