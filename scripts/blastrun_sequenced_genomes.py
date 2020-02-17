#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 10:20:27 2020

@author: meike
"""

'''
Blast run: makes bash script for the sequenced genomes and S. VT_162 (streptococcus_01814).
'''
import os
from pathlib import Path 

def blast_run_bash(db_ids, savedir):
    '''
    writes bash line to run blast for all ids in taxon_list:
    blastp -query blastquery/DB_ID.fasta  -db blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  
    -evalue 1e-5  -outfmt 6 -num_threads 8 -out blastres/DB_ID.tab
    '''
    homepath = '/home/meiker/orthomcl/'
    with open (savedir, 'w') as f:
        for id_ in db_ids:
            f.write("blastp -query " + homepath + "blastquery/"+id_+".fasta  -db " + homepath + "blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads 8 -out " + homepath + "blastres/"+id_+".tab\n")
            
path = os.getcwd()
p = Path(path)
            
db_ids =[]

with open (os.path.join(p.parents[0], 'files', '20012020streptococcus_patric_id_with_database_id.tsv')) as f:
    f.readline()
    for line in f:
        a = line.strip().split('\t')
        if a[0] == 'streptococcus_01814':
            db_ids.append(a[0])
        if len(a) >2:
            db_ids.append(a[0])


#split ids over two files
blast_run_bash(db_ids[:13], os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'surprise', 'blastrun_seq_genomes_1.sh'))

blast_run_bash(db_ids[13:], os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'surprise', 'blastrun_seq_genomes_2.sh'))