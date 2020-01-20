#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 12:07:50 2020

@author: meike
"""

'''find missing adjusted fastas'''

import os
from pathlib import Path

path = os.getcwd()
p = Path(path)


db_ids = []
missing_ids =[]
with open(os.path.join(p.parents[0],'files', 'random_files', 'adjusted_fasta_ids')) as f:
    with open(os.path.join(p.parents[0],'files', 'random_files', 'annotated_genes_ids')) as f2:
        for line in f:
            line = line.strip().split('.')
            db_ids.append(line[0])
        for line in f2:
            line = line.strip()
            if line not in db_ids:
                missing_ids.append(line)

with open (os.path.join(p.parents[0], 'scripts','bash_scripts', '20012020_missing_streptococcus_porthomcl_prep.sh'), 'w') as f:
        for id_ in missing_ids:
            f.write('orthomclAdjustFasta '+id_+' /home/meiker/git/data/prokka_annotation/'+id_+'/'+id_+'.faa 1\n')
        f.write('mv *.fasta /home/meiker/orthomcl/compliantFasta')
           