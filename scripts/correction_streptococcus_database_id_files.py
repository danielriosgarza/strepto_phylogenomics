#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 15:37:57 2020

@author: meike
"""

'''
Correction of a database_id. 

Two identical ids with same patric id, however, different genomes/species. db_id: streptococcus_04796 (instead of streptococcus_06010) with patric id 1311.146 instead of 1311.1460 (checked manually on Patric site: https://www.patricbrc.org/search/?keyword(%221311.1460%22))
'''
import os
from pathlib import Path


path = os.getcwd()
p = Path(path)
files_dir = os.path.join(p.parents[1], 'files')

#correct the streptococcus database file
with open(files_dir + '/20012020_streptococcus_database.tsv') as f:
    with open(files_dir + '/02032020_streptococcus_database.tsv', 'w') as f2:
        headers = f.readline().strip().split('\t')
        db_ind = headers.index('database_id')
        patric_ind = headers.index('genome_id')
        name_ind = headers.index('genome_name')
        f2.write('\t'.join(headers) + '\n')
        for line in f:
            a = line.strip().split('\t')
            if a[db_ind] == 'streptococcus_04796' and 'C001' in a[name_ind]:
                new_line = a
                new_line[db_ind] = 'streptococcus_06010'
                new_line[patric_ind] = '1311.1460'
                f2.write('\t'.join(new_line) + '\n')
            else:
                f2.write(line)
                
with open(files_dir + '/02032020_streptococcus_database.tsv') as f:
    headers = f.readline().strip().split('\t')
    patric_id_ind = headers.index('genome_id')
    db_id_ind = headers.index('database_id')
    species_ind = headers.index('species')
    with open(files_dir + '/02032020_streptococcus_patric_id_with_database_id.tsv', 'w') as f2:
        f2.write('\t'.join(['database_id', 'genome_id', 'species']) +'\n')
        for line in f:
            a = line.strip().split('\t')
            if a[species_ind].startswith('RB'):
                f2.write(a[db_id_ind] + '\t\t' + a[species_ind] + '\n')
            else:
                f2.write(a[db_id_ind] + '\t' + a[patric_id_ind] + '\n')