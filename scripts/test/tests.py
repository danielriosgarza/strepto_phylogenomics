#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 09:16:34 2019

@author: meike
"""
import random
from datetime import date



def split_files(db_ids):
    '''
    Splits db_list depending on length and makes lists with savdir that can be used for bash file generation.
    Tuples w/ (db_id, original_index)
    '''    
    groupsize = int(len(db_ids)/6)
    dbwithi = [(id_, i) for i, id_ in enumerate(db_ids)]
    
    id_i = [dbwithi[i :i +groupsize] for i in range(0, len(dbwithi), groupsize)]
           
    for i in range(0, len(db_ids), groupsize):
        id_i.append([db_ids(db_ids[i],db_ids.index(db_ids[i])))   
    
    [lst[i:i + n] for i in range(0, len(lst), n)]
    
    return id_i


def histogram(item_l, title=None, xlabel=None, label=None, bins=10, color='#2166ac', rotation=0):
    '''
    Makes density histogram from given list.
    '''
    plt.hist(item_l, bins, density=1, alpha=0.75, color = color, label =label)
    plt.title(title, fontsize= 14)
    plt.xlabel(xlabel, fontsize= 12)
    plt.xticks(rotation=rotation)
    plt.ylabel("Density", fontsize= 12)
    plt.tight_layout() 
    
    #plt.savefig(os.path.join(p.parents[0], 'figures', )

test = 

#%%

#run in terminal: pyhton and run following code (to check if alignments are of different length)
import os
def count_line_len(file):
    with open(file) as f:
        f.readline()
        f.readline()
        l1 = f.readline().split('       ')[1]
        for line in f:
            try:
                if len(line.split('       ')[1])!=len(l1):
                    print('different')
            except(IndexError):
              #  print (line)

for file in os.listdir('/home/meiker/phylo_tree/msa_trimmed'):
    count_line_len(file)
 
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)
print (today)

#%%
#import random
from random import shuffle
import os
from pathlib import Path
from datetime import date

def blast_Parser_bash(db_ids, savedir):
    '''
    Writes bash for BlastParser for all ids in taxon_list
    porthomclBlastParser blastres/id_.tab compliantFasta >> splitSimSeq/id_.ss.tsv
    '''
    with open (savedir, 'w') as f:
        for id_ in db_ids:
            f.write("porthomclBlastParser /home/meiker/orthomcl/blastres/"+id_+".tab /home/meiker/orthomcl/compliantFasta >> /home/meiker/orthomcl/splitSimSeq/"+id_+".ss.tsv\n")
            
def finding_best_hits(db_ids, savedir):
    '''
     paralogs are found, and an unnormalized score is assigned to them. Step 5.3 will normalize 
     this score so that it be comparable among different genomes.
    '''
    with open(savedir, 'w') as f:
        for i, id_ in enumerate(db_ids):
            f.write("porthomclPairsBestHit.py -t /home/meiker/orthomcl/taxon_list -s /home/meiker/orthomcl/splitSimSeq -b /home/meiker/orthomcl/besthit -q /home/meiker/orthomcl/paralogTemp -x "+str(i + 1)+"\n")
  
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
            f.write("porthomclPairsOrthologs.py -t /home/meiker/orthomcl/taxon_list -b /home/meiker/orthomcl/besthit -o /home/meiker/orthomcl/orthologs -x "+ str(i + 1)+"\n")

def find_paralogs(db_ids, savedir):
    '''
    Bash lines for finding paralogs. Uses split lists (db_id, index).
    '''
    with open (savedir, 'w') as f:
        for i, id_ in enumerate(db_ids):
            f.write ("porthomclPairsInParalogs.py -t /home/meiker/orthomcl/taxon_list -q /home/meiker/orthomcl/paralogTemp -o /home/meiker/orthomcl/ogenes -p /home/meiker/orthomcl/paralogs -x "+str(i + 1)+"\n")

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)

#Get ids_ of all genomes that are already blasted (located in blastres/)
dbs_done = []
for file in list(os.listdir('/home/meike/strepto_phylogenomics/files/random_files/test_updater')):    
    id_ = file.strip().split('.')[0]
    dbs_done.append(id_)
    
dbs_completed = []
for files in list(os.listdir('/home/meike/strepto_phylogenomics/files/random_files/test_updater')):
    id_ = file.strip().split('.')[0]
    dbs_done.append(id_)

blast_Parser_bash(dbs_done,os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'tests', today+'_test_blastparser.sh'))

