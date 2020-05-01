#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 10:23:00 2020

@author: meike
"""


'''
Copy chosen gff files to roary directory
'''
import os
from pathlib import Path
import statistics

path = os.getcwd()
p = Path(path)

ids = []
with open('/home/meiker/phylo_tree/roary/29042020_reduced_concat_alignments.fa') as f:
    for line in f:
        if line.startswith('>'):
            id_ = line.strip()[1::]
            ids.append(id_)
            
ids2gs = {}

with open (os.path.join(p.parents[0], 'files', '23032020_prokka_genome_data.tsv')) as f:
    f.readline()
    for line in f:
        a = line.strip().split()
        ids2gs[a[0]] = int(a[1])
        
sizes = []
for i in ids:
    sizes.append(ids2gs[i])
    
mean = statistics.mean(sizes)
dev = statistics.stdev(sizes)
lower = mean - dev
upper = mean + dev

filtered_ids = []
for i in ids:
    gs = ids2gs[i]
    if gs <= lower or gs >= upper:
        pass
    else:
        filtered_ids.append(i)

with open(os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'copy_gffs.sh'), 'w') as f:
    for i in filtered_ids:
        f.write('cp /home/meiker/git/data/prokka_annotation/' + i + '/' + i + '.gff /home/meiker/phylo_tree/roary/filtered_gff_files\n')