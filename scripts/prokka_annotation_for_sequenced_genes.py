#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 15:50:04 2020

@author: meike
"""

'''
Make prokka annotation for the 24 sequenced genomes.
'''

import os
from pathlib import Path

#prokka /home/meiker/git/genomes/streptococcus_00001.fna --outdir /home/meiker/git/prokka_annotation/streptococcus_00001 --prefix streptococcus_00001 --genus streptococcus 

path = os.getcwd()
p = Path(path)

with open(os.path.join(p.parents[0], 'files' , '17012020streptococcus_patric_id_with_database_id.tsv')) as f:
    ids = []
    species = []
    headers = f.readline().strip().split('\t')
    for line in f:
        a = line.strip().split('\t')
        if a[-1].startswith('RB'):
            ids.append(a[0])
            species.append(a[-1]+'.fasta')

with open(os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'prokka_annotation_sequenced_genomes.sh'), 'w') as f:
    for i in range(len(ids)):
        f.write('prokka /home/meiker/genomes/' +species[i]+ ' --outdir /home/meiker/git/prokka_annotation/'+ids[i]+' --prefix ' +ids[i]+' --genus streptococcus\n')