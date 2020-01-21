#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 12:10:37 2020

@author: meike
"""

'''
Run first the porthoMCL_prep for all species and follow than the steps in the PorthoPrep Steps document in the terminal,
before running the rest.


Split funciton not ready!
'''
import random
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
    homepath = '/home/meiker/orthomcl/'
    with open (savedir, 'w') as f:
        for id_ in db_ids:
            f.write("blastp -query " + homepath + "blastquery/"+id_[0]+".fasta  -db " + homepath + "blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads 8 -out " + homepath + "blastres/"+id_[0]+".tab\n")
    
def blast_Parser_bash(db_ids, savedir):
    '''
    Writes bash for BlastParser for all ids in taxon_list
    porthomclBlastParser blastres/id_.tab compliantFasta >> splitSimSeq/id_.ss.tsv
    '''
    with open (savedir, 'w') as f:
        for id_ in db_ids:
            f.write("porthomclBlastParser blastres/"+id_+".tab compliantFasta >> splitSimSeq/"+id_+".ss.tsv\n")
            
def finding_best_hits(db_ids, savedir):
    '''
     paralogs are found, and an unnormalized score is assigned to them. Step 5.3 will normalize 
     this score so that it be comparable among different genomes.
    '''
    with open(savedir, 'w') as f:
        for i, id_ in enumerate(db_ids):
            f.write("porthomclPairsBestHit.py -t taxon_list -s splitSimSeq -b besthit -q paralogTemp -x "+str(i + 1)+"\n")
  
def split_files(db_ids):
    '''
    Splits db_list depending on length and makes lists with savdir that can be used for bash file generation.
    Tuples w/ (db_id, original_index)
    '''    
    #groupsize = int(len(db_ids)/16)
    groupsize = 761
    db_index = [(id_, i + 1) for i, id_ in enumerate(db_ids)]
    
    id_i = [db_index[i :i + groupsize] for i in range(0, len(db_index), groupsize)]
    
    return id_i

def find_orthologs(db_ids, savedir):
    '''
    Output of bash line is all the ortholog genes.
    '''
    with open (savedir, 'w') as f:
        for i, id_ in enumerate(db_ids):
            f.write("porthomclPairsOrthologs.py -t taxon_list -b besthit -o orthologs -x "+ str(i + 1)+"\n")

def find_paralogs(db_ids, savedir):
    '''
    Bash lines for finding paralogs. Uses split lists (db_id, index).
    '''
    with open (savedir, 'w') as f:
        for i, id_ in enumerate(db_ids):
            f.write ("porthomclPairsInParalogs.py -t taxon_list -q paralogTemp -o ogenes -p paralogs -x "+str(i + 1)+"\n")

    
path = os.getcwd()
p = Path(path)


#blast_run_bash(db_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20200114_floricoccus_blastrun.sh'))
#blast_Parser_bash(db_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20200114_floricoccus_blastparser.sh'))
#finding_best_hits(db_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20200114_floricoccus_find_best_hits.sh'))

#get ids
strepto_ids = get_ids(os.path.join(p.parents[0], 'files', '20012020streptococcus_patric_id_with_database_id.tsv'))
flori_ids = get_ids(os.path.join(p.parents[0], 'files', 'floricoccus_patric_id_with_database_id.tsv'))
lacto_ids = get_ids(os.path.join(p.parents[0], 'files', 'lactococcus_patric_id_with_database_id.tsv'))

#adjust fasta files
porthoMCL_prep(strepto_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20012020_streptococcus_porthomcl_prep.sh'))
porthoMCL_prep(flori_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20012020_floricoccus_porthomcl_prep.sh'))
porthoMCL_prep(lacto_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', '20012020_lactococcus_porthomcl_prep.sh'))

#get all ids in a single list and split it for blast run
db_ids = get_taxon_list(os.path.join(p.parents[0], 'files', 'taxon_list'))
splitted_ids = split_files(db_ids)

for i, l_ids in enumerate(splitted_ids):
    blast_run_bash(l_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts' , 'porthomcl', '200120_blastrun'+str(i)+'.sh'))


blast_Parser_bash(db_ids,os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', '210120_blastparser.sh'))
finding_best_hits(db_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', '210120_find_best_hits.sh'))
find_orthologs(db_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', '210120_orthologs.sh'))
find_paralogs(db_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', '210120_paralogs.sh'))


# test_ids = get_taxon_list(os.path.join(p.parents[0], "files", 'porthomcl', 'taxon_list'))
# ids_split = split_files(test_ids)

# for i, l_ids in enumerate(ids_split):
#     blast_run_bash(l_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts' , 'porthomcl', 'testset_blastrun'+str(i)+'.sh'))

# blast_Parser_bash(test_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl','testset_blastparser.sh'))
# finding_best_hits(test_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl','testset_besthits.sh'))
# find_orthologs(test_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'testset_orthologs.sh'))
# find_paralogs(test_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'testset_paralogs.sh'))    