#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 15:01:14 2020

@author: meike
"""


'''Run missing best hits'''

import os
from pathlib import Path
from datetime import date
   
path = os.getcwd()
p = Path(path)

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)

ids_done = []

for file in list(os.listdir('/home/meiker/orthomcl/besthit')):    
    id_ = file.strip().split('.')[0]
    ids_done.append(id_)

missing = []
for file in list(os.listdir('/home/meiker/orthomcl/splitSimSeq/')):
    id_ = file.strip().split('.')[0]
    if id_ not in ids_done:
        missing.append(id_)
        
taxon_list = []
with open(os.path.join(p.parents[0], 'files', 'taxon_list')) as f:
    for line in f:
        taxon_list.append(line.strip())

indexes = []        
for i, taxon in enumerate(taxon_list):
    if taxon in missing:
        indexes.append(i+1)

print(len(missing))
    
size = round(len(indexes)/3)

for i in range(1,4):
    with open(os.path.join(p.parents[0], 'scripts', 'bash_scripts', 'porthomcl', 'besthit', + today + '_missing_best_hits' + str(i) + '.sh'), 'w') as f:
        for j in indexes:
            f.write("porthomclPairsBestHit.py -t /home/meiker/orthomcl/taxon_list -s /home/meiker/orthomcl/splitSimSeq -b /home/meiker/orthomcl/besthit -q /home/meiker/orthomcl/paralogTemp -x "+str(j)+"\n")