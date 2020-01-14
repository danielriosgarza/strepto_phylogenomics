#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 12:10:37 2020

@author: meike
"""

'''
Run first the porthoMCL-prep for all species and follow than the steps in the PorthoPrep document in the terminal.

PorthoMCL preparation, when prep runned: 
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

def get_taxon_list(taxon_list):
    '''
    After creation of taxon_list in terminal, run function to get all ids.
    '''
    db_ids = []
    with open(taxon_list) as f:
        for line in f:
            line = line.strip()
            db_ids.append(line)
    return db_ids

def blast_run_bash(db_ids, savedir):
    '''
    writes bash line to run blast for all ids in taxon_list:
    blastp -query blastquery/DB_ID.fasta  -db blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  
    -evalue 1e-5  -outfmt 6 -num_threads 8 -out blastres/DB_ID.tab
    '''
    with open (savedir, 'w') as f:
        for id_ in db_ids:
            f.write("blastp -query blastquery/"+id_+".fasta  -db blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads 8 -out blastres/"+id_+".tab\n")
    
def blast_Parser_bash(db_ids, savedir):
    '''
    Writes bash for BlastParser for all ids in taxon_list
    porthomclBlastParser blastres/id_.tab compliantFasta >> splitSimSeq/id_.ss.tsv
    '''
    with open (savedir, 'w') as f:
        for id_ in db_ids:
            f.write("porthomclBlastParser blastres/"+id_+".tab compliantFasta >> splitSimSeq/"+id_+".ss.tsv\n")
            


path = os.getcwd()
p = Path(path)

flori_ids = get_ids(os.path.join(p.parents[0], 'files', 'floricoccus_patric_id_with_database_id.tsv'))
sav_test = os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20200113_floricoccus_PorthoMCL_prep.sh')

porthoMCL_prep(flori_ids, sav_test)

db_ids = get_taxon_list(os.path.join(p.parents[0], "files", 'porthomcl', 'taxon_list'))

blast_run_bash(db_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20200114_floricoccus_blastrun.sh'))
blast_Parser_bash(db_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20200114_floricoccus_blastparser.sh'))