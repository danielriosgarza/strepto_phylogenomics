#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 11:05:15 2020

@author: meike
"""

'''
Swiftortho test. Random files from 10 streptococcus genomes.
'''
import random
import os
from pathlib import Path

path = os.getcwd()
p = Path(path)


db_ids = []
with open(os.path.join(p.parents[0], 'files', 'taxon_list')) as f:
    for line in f:
        line = line.strip()
        db_ids.append(line)
db_ids

#Get testset with random ids
test = random.sample(db_ids, 10)

with open(os.path.join(p.parents[0],'scripts', 'bash_scripts', 'blat', '280120_pblat_run.sh'), 'w') as f:
    for id_ in test:
        f.write('/home/meiker/software/icebert-pblat-652d3b3/pblat -prot -threads=8 -out=blast /home/meiker/orthomcl/filteredFasta/goodProteins.fasta /home/meiker/git/data/prokka_annotation/'+id_+'/'+id_+'.faa /home/meiker/pblat/'+id_+'.blast\n\n')