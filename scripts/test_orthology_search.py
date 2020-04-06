#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 10:34:21 2020

@author: meike
"""
'''
Testing of the orthology search. Make taxon list that only includes the already parsed ids and make the following steps depending on that new taxon list
'''

import os
from pathlib import Path
from datetime import date

def finding_best_hits(indexes_ids, savedir):
    '''
     paralogs are found, and an unnormalized score is assigned to them. Step 5.3 will normalize 
     this score so that it be comparable among different genomes.
    '''
    with open(savedir, 'w') as f:
        for i in indexes_ids:
            f.write("porthomclPairsBestHit.py -t /home/meiker/orthomcl/taxon_list2 -s /home/meiker/orthomcl/splitSimSeq -b /home/meiker/orthomcl/besthit2 -q /home/meiker/orthomcl/paralogTemp2 -x "+str(i)+" -l /home/meiker/orthomcl/besthit2\n")

def split_files(i_list, nsplits = 4):
    '''
    Splits list depending on length.
    '''    
    groupsize = round(len(i_list)/nsplits)
    
    inds = [i_list[i:i + groupsize] for i in range(0, len(i_list), groupsize)]
    return inds

#get all ids (genomes) that were already blasted and parsed
path = os.getcwd()
p = Path(path)

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)

dbs_parsed = []

for file in list(os.listdir('/home/meiker/orthomcl/splitSimSeq/')):
    id_ = file.strip().split('.')[0]
    dbs_parsed.append(id_)
    
dbs_parsed = sorted(dbs_parsed)

indexes = []
with open ('/home/meiker/orthomcl/taxon_list2' , 'w') as f:
    for i, id_ in enumerate(dbs_parsed):
        i += 1
        f.write(id_ + '\n')
        indexes.append(i)
        
inds_l = split_files(indexes)

for i, inds in enumerate(inds_l):
    finding_best_hits(inds, os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'besthit', today + '_test_find_best_hits' + str(i) + '.sh'))       


line_orthosearch = 'porthomclPairsOrthologs.py -l -t /home/meiker/orthomcl/taxon_list2 -b /home/meiker/orthomcl/besthit2 -o /home/meiker/orthomcl/orthologs2 -x '