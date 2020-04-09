#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 07:58:43 2020

@author: meike
"""


'''
Testing of the Prothomcl protocol with a small subset of 50 random genomes.
'''
from datetime import date
import os
from pathlib import Path
from random import sample

def blast_Parser_bash(ids, savedir):
    '''
    Writes bash for BlastParser for all ids in taxon_list
    porthomclBlastParser blastres/id_.tab compliantFasta >> splitSimSeq/id_.ss.tsv
    '''
    with open (savedir, 'w') as f:
        for id_ in ids:
            f.write("porthomclBlastParser /home/meiker/tests/orthomcl/blastres/" + id_[0] + ".tab /home/meiker/tests/orthomcl/compliantFasta >> /home/meiker/tests/orthomcl/splitSimSeq/" + id_[0] + ".ss.tsv\n")
            
def blast_run_bash(ids, savedir):
    '''
    writes bash line to run blast for all ids in taxon_list:
    blastp -query blastquery/DB_ID.fasta  -db blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  
    -evalue 1e-5  -outfmt 6 -num_threads 8 -out blastres/DB_ID.tab
    '''
    homepath = '/home/meiker/tests/orthomcl/'
    with open (savedir, 'w') as f:
        for id_ in ids:
            f.write("blastp -query " + homepath + "blastquery/" + id_[0] + ".fasta  -db " + homepath + "blastdb/goodProteins.fasta  -seg yes  -dbsize 100000000  -evalue 1e-5  -outfmt 6 -num_threads 8 -out " + homepath + "blastres/" + id_[0] + ".tab\n")
            
def finding_best_hits(taxons, savedir):
    '''
     paralogs are found, and an unnormalized score is assigned to them. Step 5.3 will normalize 
     this score so that it be comparable among different genomes.
    '''
    with open(savedir, 'w') as f:
        for id_ in enumerate(taxons):
            f.write('porthomclPairsBestHit.py -t /home/meiker/tests/orthomcl/taxon_list -s /home/meiker/tests/orthomcl/splitSimSeq -b /home/meiker/tests/orthomcl/besthit -q /home/meiker/tests/orthomcl/paralogTemp -x ' + str(id_[1]) + ' -l /home/meiker/tests/orthomcl/logs/' + today + '_logfile_besthits.txt\n')
            
def find_orthologs(taxons, savedir):
    '''
    Output of bash line is all the ortholog genes.
    '''
    with open (savedir, 'w') as f:
        for id_ in taxons:
            f.write("porthomclPairsOrthologs.py -t /home/meiker/tests/orthomcl/taxon_list -b /home/meiker/tests/orthomcl/besthit -o /home/meiker/tests/orthomcl/orthologs -x " + str(id_[1]) + " -l /home/meiker/tests/orthomcl/logs/" + today + "_logfile_orthologs.txt\n")

def porthoMCL_prep(ids, savedir):
    '''
    Writes bash lines for PorthoMCL preparation. 
    '''
    with open (savedir, 'w') as f:
        for id_ in ids:
            f.write('orthomclAdjustFasta '+id_+' /home/meiker/git/data/prokka_annotation/' + id_ + '/' + id_ + '.faa 1\n')
        f.write('mv *.fasta /home/meiker/tests/orthomcl/compliantFasta')          
        
def split_files(ids, nsplits = 5):
    '''
    Splits db_list depending on length and makes lists with savdir that can be used for bash file generation.
    Tuples w/ (db_id, original_index)
    '''    
    groupsize = round(int(len(ids)/nsplits))
    
    db_index = [(id_, i + 1) for i, id_ in enumerate(ids)]
    
    id_i = [db_index[i :i + groupsize] for i in range(0, len(db_index), groupsize)]
    
    return id_i
               
path = os.getcwd()
p = Path(path)

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)


# taxons = []
# #look in the taxon list for random ids that are used for the porthomcl
# with open (os.path.join(p.parents[0], 'files', 'taxon_list')) as f:
#     for line in f:
#         taxons.append(line.strip())
        
# test_set = sample(taxons, 50)


# porthoMCL_prep(test_set, os.path.join(p.parents[0], 'scripts', 'bash_scripts', today + 'testset_porthomcl_prep.sh'))

# splitted_ids = split_files(test_set)

# for i, l_ids in enumerate(splitted_ids):
#     blast_run_bash(l_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts' , 'porthomcl', today + 'testset_blastrun' + str(i)+ '.sh'))

taxons = []
indexes = []
with open('/home/meiker/tests/orthomcl/taxon_list') as f:
    i = 0
    for line in f:
        i += 1
        taxons.append(line.strip())
        indexes.append(i)

splitted_ids = split_files(taxons)
    
for i, l_ids in enumerate(splitted_ids):
    i += 1
    blast_Parser_bash(l_ids, os.path.join(p.parents[0], 'scripts', 'bash_scripts' , 'porthomcl', 'blastparser', today + '_blastparser' + str(i) + '.sh'))
    

for i, l_inds in enumerate(splitted_ids):
    i += 1
    finding_best_hits(l_inds, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'besthit', today + '_find_best_hits' + str(i) +'.sh'))
    
for i, l_inds in enumerate(splitted_ids):
    i += 1
    find_orthologs(l_inds, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'orthologs', today + '_orthologs' + str(i) +'.sh'))

