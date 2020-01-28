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

with open(os.path.join(p.parents[0],'scripts', 'bash_scripts', 'blat', '280120_blat_run.sh'), 'w') as f:
    for id_ in test:
        f.write('/home/meiker/software/blat -prot /home/meiker/orthomcl/filteredFasta/goodProteins.fasta /home/meiker/git/data/prokka_annotation/' +id_+'/'+id_+'.faa /home/meiker/blat/'+id_+'.blast -out=blast -threads=8\n\n')