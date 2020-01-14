#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 12:10:37 2020

@author: meike
"""

'''
PorthoMCL preparation, when bash script runned: 
$ls -1 ~/orthomcl/compliantFasta/ | sed -e 's/\..*$//'  > taxon_list
(in orthomcl dir to create taxon_list)
'''

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
            f.write('orthomclAdjustFasta '+id_+' /home/meiker/git/data/prokka_annotation/'+id_+'/'+id_+'.faa 1\n')
        f.write('mv *.fasta /home/meiker/orthomcl/compliantFasta')

def blast_run_bash(db_ids, savedir):
    '''
    writes bash line to run blast for each db_id:
    blastp -query blastquery/DB_ID.fasta  -db blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  
    -evalue 1e-5  -outfmt 6 -num_threads 8 -out blastres/DB_ID.tab
    '''
    with open (savedir, 'w') as f:
        for id_ in db_ids:
            f.write("blastp -query blastquery/"+id_+".fasta  -db blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads 8 -out blastres/"+id_+".tab\n")
    


path = os.getcwd()
p = Path(path)

flori_ids = get_ids(os.path.join(p.parents[0], 'files', 'floricoccus_patric_id_with_database_id.tsv'))
sav_test = os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20200113_floricoccus_PorthoMCL_prep.sh')

porthoMCL_prep(flori_ids, sav_test)
blast_run_bash(flori_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20200114_floricoccus_blastrun.sh'))
