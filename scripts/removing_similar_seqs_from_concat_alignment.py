#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 13:49:04 2020

@author: meike
"""

'''
Look through MSA concat file for sequences with 95% similarity and remove the "duplicates"
'''

import os
from pathlib import Path
import itertools
from fuzzywuzzy import fuzz


path = os.getcwd()
p = Path(path)

concatfile = '/home/meiker/phylo_tree/iqtree/concat_alignments'
#concat_file = '/home/meike/tests/Files/concat_alignments'
#testfile = '/home/meike/tests/Files/prottest_testset1.fa'

seqs = {}
with open(concatfile) as f:
    name = ''
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            name = line[1:]
            seqs[name] = ''
        else:
            seqs[name] += line

# generate a map of unique sequence to identical sequences
similar = []

names = []
sequences = []
for name, seq in seqs.items():
    names.append(name)
    sequences.append(seq)

#ratios = []
    
#iterate over the sequences and compare each sequence with each other
for a, b in itertools.combinations(sequences, 2):
    Ratio = fuzz.ratio(a,b)
    #print(count, ' Ratio determined')
    
    #if the sequences are at least 95% identical add the ids (names) to the list 'similar'
    if Ratio >= 95:
        #ratios.append(Ratio)
        a_ind = sequences.index(a)
        b_ind = sequences.index(b)
        similar.append((names[a_ind], names[b_ind]))

#make dict to store ids (first appeared) as keys and all ids with similar seqs as values
ids2dupids = {}

left_over_ids = []
for i in similar:
    #check if id is already a key or values (list of duplicate ids)
    if i[0] not in ids2dupids.keys() and i[0] not in [x for v in ids2dupids.values() for x in v]:
        ids2dupids[i[0]] = [i[1]]
        left_over_ids.append(i[0])
    else:
        #loop through dict and look at keys and duplicates
        for key, duplicate in ids2dupids.items():
            #if the first item of the similar sequences is a key, than add the second to the values as duplicate
            if i[0] == key:
                ids2dupids[key].append(i[1])
print("Number of left over IDs: ", len(left_over_ids),"\nleft over IDs: ", left_over_ids)
                

unique_seq_ids = ids2dupids.keys()
              
with open('/home/meiker/phylo_tree/iqtree/reduced_alignments/reduced_concat_alignments.fa', 'w') as f:
    for k in unique_seq_ids:
        samp = (k, seqs[k])
        f.write(">{}\n{}\n".format(*samp))
            

            
