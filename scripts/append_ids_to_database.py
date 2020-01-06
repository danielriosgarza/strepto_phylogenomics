#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 09:19:34 2020

@author: daniel
"""

import os
from pathlib import Path

path = os.getcwd()
p = Path(path)

def append_database_ids(downloads_file, database_file, new_database_file):

    names_db ={}
    c=0
    with open(downloads_file) as f:
        
        for line in f:
            a=line.strip().split(' ')
            
            names_db[a[1].split('/')[-2]] =a[3].split('/')[-1].replace('.fna', '')
    
    with open(database_file) as f:
        a = f.readline().strip().split('\t')
        gidind = a[1::].index('genome.genome_id')
        
        with open(new_database_file, 'w') as f2:
            
            a = ['database_id'] + a[1::]
            f2.write('\t'.join(a)+'\n')
            for line in f:
                a=line.strip().split('\t')
                while a[0] not in names_db:
                    a[0] = a[0]+'0'
                    a[gidind]=a[0]
                gid = names_db[a[0]]
                a = [gid] + a[1::]
                f2.write('\t'.join(a)+'\n')
        

append_database_ids(os.path.join(p,'bash_scripts',  'get_lactococcus_genomes_patric.sh'), os.path.join(p.parents[0], 'files', 'lactococcus_genomes_quality.tsv'),os.path.join(p.parents[0], 'files', '06012020_lactococcus_genomes_quality.tsv'))
append_database_ids(os.path.join(p,'bash_scripts',  'get_floricoccus_genomes_patric.sh'), os.path.join(p.parents[0], 'files', 'floricoccus_genomes_quality.tsv'),os.path.join(p.parents[0], 'files', '06012020_floricoccus_genomes_quality.tsv'))
append_database_ids(os.path.join(p, 'bash_scripts', 'get_streptococcus_genomes_patric.sh'), os.path.join(p.parents[0], 'files', 'streptococcus_genomes_quality.tsv'),os.path.join(p.parents[0], 'files', '06012020_streptococcus_genomes_quality.tsv'))


