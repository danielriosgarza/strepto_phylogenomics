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
test_swiftortho = random.sample(db_ids, 10)

with open(os.path.join(p.parents[0],'scripts', 'bash_scripts', 'swiftortho', 'homologous_search_test.sh'), 'w') as f:
    for id_ in test_swiftortho:
        f.write('python /home/meiker/software/SwiftOrtho/bin/find_hit.py -p blastp -i /home/meiker/git/data/prokka_annotation/'+id_+'/'+id_+'.faa -d /home/meiker/orthomcl/filteredFasta/goodProteins.fasta -o /home/meiker/swiftortho/'+id_+'.faa -e 1e-5 -s 111111 -a 8\n\n')