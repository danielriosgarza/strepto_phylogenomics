#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 11:17:46 2020

@author: meike
"""

'''Getting genome sizes of sequenced genomes and add them to a table of information of the sequenced genomes'''

import os
from pathlib import Path
    
path = os.getcwd()
p = Path(path)


first_id = 'streptococcus_11962'

seq_ids = []
        
with open (os.path.join(p.parents[0], 'files', '20012020streptococcus_patric_id_with_database_id.tsv')) as f:
    seen = False
    for line in f:
        a = line.strip().split('\t')
        if a[0] == first_id:
            seen = True
        if seen:
            seq_ids.append(a[0])

#Look in the log files from prokka annotation to get genome sizes            
gs = {}
for folder in os.scandir('/home/meiker/git/data/prokka_annotation'):
    _id = str(folder).split("'")[1]
    if _id in seq_ids:
        with open(folder.path + '/' + _id + '.log') as f:
            for line in f:
                if line.startswith('Contigs'):
                    a = line.strip().split(' ')
                    gs[_id] = a[2]
                    
with open(os.path.join(p.parents[0], 'files', '23032020_sequenced_genomes_database.tsv'), 'w') as f:
    f.write('database_id\tgenome_length\n')
    for i in seq_ids:
        f.write(i + '\t' + gs[i] + '\n')



      