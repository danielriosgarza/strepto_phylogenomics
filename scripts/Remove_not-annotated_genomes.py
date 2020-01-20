#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 10:09:57 2020

@author: meike
"""

'''Remove #11895 and #11939 from databases, because they lack .fna files --> no prokka annotation'''

# 17012020streptococcus_patric_id_with_database_id.tsv
# 17012020_streptococcus_database.tsv

import os
from pathlib import Path

path = os.getcwd()
p = Path(path)

db_table = os.path.join(p.parents[0], 'files', '17012020streptococcus_patric_id_with_database_id.tsv')
meta_table = os.path.join(p.parents[0], 'files', '17012020_streptococcus_database.tsv')

with open (db_table) as f:
    with open(os.path.join(p.parents[0], 'files', '20012020streptococcus_patric_id_with_database_id.tsv'), 'w') as f2:
        headers = f.readline()
        f2.write(headers)
        for line in f:
            a = line.strip().split('\t')
            if a[0] != "streptococcus_11895" and a[0] != 'streptococcus_11939':
                f2.write('\t'.join(a)+'\n')

with open(meta_table) as f:
    with open(os.path.join(p.parents[0], 'files', '20012020_streptococcus_database.tsv'), 'w') as f2:
        headers = f.readline()
        f2.write(headers)
        for line in f:
            a = line.strip().split('\t')
            if a[0] != "streptococcus_11895" and a[0] != 'streptococcus_11939':
                f2.write('\t'.join(a)+'\n')