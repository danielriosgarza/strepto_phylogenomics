#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 12:10:37 2020

@author: meike
"""

'''
PorthoMCL preparation
'''
#orthomclAdjustFasta NC_000913 sample/0.input_faa/NC_000913.faa 4
#
#orthomclAdjustFasta label/input_fasta/number_identifiying_field_containing_fasta_header

#orthomclAdjustFasta database_id/db_id.faa/last_column.tsv ?

import os
from pathlib import Path

def get_ids(file):
    '''
    Gets ids from database-patric_id file and returns a list with all ids
    '''
    with open (file) as f:
        db_ids = []
        for line in f:
            a = line.strip().split('\t')
            if a[0] != 'database_id':
                db_ids.append(a[0])
        db_ids.sort()
    return db_ids

def porthoMCL_prep(db_ids, savedir):
    '''
    Writes bash lines for PorthoMCL preparation. 
    '''
    with open (savedir, 'w') as f:
        for id_ in db_ids:
            f.write('orthomclAdjustFasta '+id_+' /home/meiker/git/data/prokka_annotation/'+id_+'.faa 1\n')
        f.write('mv *.fasta /home/meiker/orthomcl/compliantFasta')
 

           
def taxon_list(db_ids, savedir):
    '''
    Makes taxon_list for 
    '''
    with open (savedir, 'w') as f:        
        for id_ in db_ids:
            f.write(id_+"\n")


path = os.getcwd()
p = Path(path)

flori_ids = get_ids(os.path.join(p.parents[0], 'files', 'floricoccus_patric_id_with_database_id.tsv'))
sav_test = os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20200113_floricoccus_PorthoMCL_prep.sh')

porthoMCL_prep(flori_ids, sav_test)
taxon_list(flori_ids, os.path.join(p.parents[0], 'files', 'taxon_lists', '20200113_floricoccus_taxon_list'))