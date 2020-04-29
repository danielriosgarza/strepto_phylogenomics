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

path = os.getcwd()
p = Path(path)

ids = []
with open('/home/meiker/phylo_tree/roary/29042020_reduced_concat_alignments.fa') as f:
    for line in f:
        if line.startswith('>'):
            id_ = line.strip()[1::]
            ids.append(id_)

with open(os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'move_gffs.sh'), 'w') as f:
    for i in ids:
        f.write('cp /home/git/data/prokka_annotation/' + i + '/' + i + '.gff /home/meiker/roary/gff_files\n')