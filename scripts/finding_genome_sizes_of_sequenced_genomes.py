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
species = []       
with open (os.path.join(p.parents[0], 'files', '20012020streptococcus_patric_id_with_database_id.tsv')) as f:
    seen = False
    for line in f:
        a = line.strip().split('\t')
        if a[0] == first_id:
            seen = True
        if seen:
            seq_ids.append(a[0])
            species.append(a[2])

#Look in the log files from prokka annotation to get genome sizes            
gs = {}

for folder in os.scandir('/home/meiker/git/data/prokka_annotation'):
    _id = str(folder).split("'")[1]
    if _id in seq_ids:
        with open(folder.path + '/' + _id + '.txt') as f:
            for line in f:
                a = line.strip().split(' ')
                if line.startswith('bases'):
                    gs[_id] = a[-1]
                
with open(os.path.join(p.parents[0], 'files', '23032020_sequenced_genomes_database2.tsv'), 'w') as f:
    f.write('database_id\tspecies\tgenome_length\tcollection_year\thost_name\tisolation_country\tisolation_source\tsequencing_platform\n')
    for i in seq_ids:
        f.write(i + '\t' + species[seq_ids.index(i)] + '\t' + gs[i] + '\n')

all_ids = []
gs = {}
CDS = {}
rRNA = {}
repeat_region = {}
tRNA = {}
for folder in os.scandir('/home/meiker/git/data/prokka_annotation'):
    _id = str(folder).split("'")[1]
    all_ids.append(_id)
    with open(folder.path + '/' + _id + '.txt') as f:
        for line in f:
            a = line.strip().split(' ')
            if line.startswith('bases'):
                gs[_id] = a[-1]
            elif line.startswith('CDS'):
                CDS[_id] = a[-1]
            elif line.startswith('rRNA'):
                rRNA[_id] = a[-1]
            elif line.startswith('repeat_region'):
                repeat_region[_id] = a[-1]
            elif line.startswith('tRNA'):
                tRNA[_id] = a[-1]

print(all_ids) 

with open(os.path.join(p.parents[0], 'files', '23032020_prokka_genome_data.tsv'), 'w') as f:
    f.write('database_id\tgenome_size\tCDS\trRNA\trepeat_region\ttRNA\n')
    for i in all_ids:
        f.write(i + '\t' + gs[i] + '\t' + CDS[i] + '\t' + rRNA[i] + '\t' + repeat_region[i] + '\t' + tRNA[i] + '\n')