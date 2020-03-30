#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 11:55:25 2020

@author: meike
"""


'''
Prepare Binary table of Pan-genome. Filter out the paralogs from the orthologs file.
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

        
#make prep table: all orthologs (prot ids) in a list with all genomes with 0 and 1 following
   
with open (os.path.join(p.parents[0], 'files', 'binary_table', 'all.ort.group')) as f:
    with open(os.path.join(p.parents[0], 'files', 'binary_table', today + '_binary_table_prep1.tsv'), 'w') as f2:
        f2.write('pan_genome\torthologs\t')
        for db_id in taxon_list:
            if db_id == taxon_list[-1]:
                f2.write(db_id + '\n')
            else:
                f2.write(db_id + '\t')
        for line in f:
            a = line.strip().split('\t')
            ids = []
            orthos = []
            for pair in a:
                id_ = pair.split('|')[0]
                ids.append(id_)
                ortho = pair.split('|')[1]
                orthos.append(ortho)
            group = [str(0)]*len(taxon_list)    
            if len(set(ids)) > 5:
                f2.write('\t')    
                #mark found ids with an 1 in the table, rest is 0
                f2.write(','.join(orthos) + '\t')
                for id_ in set(ids):
                    group[taxon_list.index(id_)] = str(1)
                f2.write('\t'.join(group) + '\n')
                
            
                
            
        