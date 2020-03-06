#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 12:18:43 2020

@author: meike
"""

'''Making Table with the patric_id and our databse id'''

import os
from pathlib import Path
from datetime import date


def patric_and_db_ids(file, genus, Date):
    with open(file) as f:
        headers = f.readline().strip().split('\t')
        headers_inds = {i : name for name, i in enumerate(headers)} 
        patric_id_ind = headers_inds['genome.genome_id']
        db_id_ind = headers_inds['database_id']
        with open(os.path.join(p.parents[0], 'files', Date + '_' + genus +  '_patric_id_with_database_id.tsv'), 'w') as f2:
            f2.write('\t'.join(['database_id', 'genome_id']) +'\n')
            for line in f:
                a = line.strip().split('\t')
                patric_id = a[patric_id_ind]
                db_id = a[db_id_ind]
                
                f2.write('\t'.join([db_id, patric_id]) + '\n')

path = os.getcwd()
p = Path(path)

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)
                
lactoinput = os.path.join(p.parents[0], 'files', '06012020_lactococcus_genomes_quality.tsv')
floriinput = os.path.join(p.parents[0], 'files', '06012020_floricoccus_genomes_quality.tsv')
streptoinput = os.path.join(p.parents[0], 'files', 'streptococcus_genomes_quality.tsv')


#patric_and_db_ids(lactoinput, 'lactococcus')
#patric_and_db_ids(floriinput, 'floricoccus')
patric_and_db_ids(streptoinput, 'streptococcus', today)