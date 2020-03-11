#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 13:49:04 2020

@author: meike
"""

'''
Look through MSA concat file for sequences with 95% similarity and remove the "duplicates"
'''

import scipy.spatial as sps
import os
from pathlib import Path
import numpy as np
from datetime import date


def replace_letter_for_number(replace_d, d_arrays):
    '''
    Needs a dictionary of letters mapping to numbers to replace and a dict with ids as keys and arrays (all letter apart) as values.
    Replaces letter for numbers. These can be compared to get the distance (similarity) of sequences.
    '''
    
    for k in d_arrays:
        #id maps to array and letter is replaced by the digit given in the replace dict.
        d_arrays[k] = np.array([replace_d[d_arrays[k][i]] for i in range(len(d_arrays[k]))])
    


def compare(idx, list_idx, v_dict, threshold = 0.05):
    '''
    Needs:
        >An id where the sequence should be compared to the rest: idx
        >A list of ids where the idx is compared to: list_idx
        >A dict having ids as keys and array with numbers instead of letters as values to determine distances (similarity) of sequences: v_dict
        >A threshold that determines when seqs are similar enough: threshold for DISTANCE (default 0.05, meaning similarity at least 95%)
        
    Returns:
        list of ids that are similar/equal to the id given (similarity depends on given threshold)
        
    Uses hamming metric to determine distances.
    '''
    v1 = [v_dict[idx]]
    v2 =np.array([v_dict[i] for i in list_idx])
    distance = sps.distance.cdist(v1,v2, metric='hamming')[0]
    
    return [list_idx[i[0]] for i in enumerate(distance) if i[1]<threshold]



path = os.getcwd()
p = Path(path)

#get the date to keep track of the scripts (added to scriptname)
today = date.today().strftime("%d/%m/%Y")
today = today.split('/')
today = ''.join(today)

concatfile = '/home/meiker/phylo_tree/iqtree/concat_alignments'
#concat_file = '/home/meike/tests/Files/concat_alignments'
#testfile = '/home/meike/tests/Files/prottest_testset1.fa'
#testfile = '/home/meike/tests/Files/prottest_testset_larger.fa'


#make dict with ids as keys and seqs as values
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

#make dict to replace letter by numbers
replace = {}
alphabet = 'abcdefghijklmnopqrstuvwxyz-'
for i, letter in enumerate(alphabet):
    replace[letter.upper()] = i+1

#put seqs in an array in order to keep number at right position (as z - 26)
seqs_array = {i:np.array(list(seqs[i])) for i in seqs}
#replace letter with numbers in order to look for similarity based on hamming distance
replace_letter_for_number(replace, seqs_array)
    
#make two lists with all ids (keys): one to iterate over and the other to update (remove already found duplicates)
k1 = list(seqs.keys())
k2 = k1[:]

#list storing all unique ids
survived = []

for i in k2:
    
    if i in k1:
        #update already the ids list that the id is not compared to itself and save unique id in survived list
        survived.append(i)
        k1.remove(i)
        #stop if k1 has no ids left
        if len(k1)==0:
            break
        #look for similar seqs and remove them from k1
        equal = compare(i, k1, seqs_array)
        for i in equal:
            k1.remove(i)

#add the in nijmegen sequenced genomes and Streptococcus_sp_VT_162 to the reduced file
with open(os.path.join(p.parents[0], 'files', '03032020_streptococcus_database_final.tsv')) as f:
    headers = f.readline().strip().split('\t')
    species_ind = headers.index('species')
    for line in f:
        a = line.strip().split('\t')
        species = a[species_ind]
        number = int(a.split('_')[1])
        if a[0] == 'streptococcus_01814' or number > 11961:
            survived.append(a[0])
        

#write all unique seqs into a concatfile that can be used to build the tree
with open('/home/meiker/phylo_tree/iqtree/reduced_alignments/rooted/' + today + '_reduced_concat_alignments.fa', 'w') as f:
    for k in survived:
        samp = (k, seqs[k])
        f.write(">{}\n{}\n".format(*samp))
        
#specify the outgroup and save it in the file 'outgroups'  
outgroup = []
strepto_ids = []
for id_ in survived:
    if 'strepto' in id_:
        strepto_ids.append(id_)
    else:
        outgroup.append(id_)
        
#file structure: ((streptococcus_00001,streptococcus_00002,...), floricoccus_00001,...);            
with open('/home/meiker/phylo_tree/iqtree/reduced_alignments/rooted/' + today + '_reduced_outgroups.txt', 'w') as f:
    f.write('((')
    for id_ in strepto_ids:
        if id_ == strepto_ids[-1]:
            f.write(id_)
        else:
            f.write(id_ + ',')
    f.write('),')
    for _id in outgroup:
        if _id == outgroup[-1]:
            f.write(_id)
        else:
            f.write(_id + ',')
    f.write((');'))


#%% runcell 2
            
###first try for removing seqs (takes very long)
            
# import os
# from pathlib import Path
# import itertools
# from fuzzywuzzy import fuzz


# path = os.getcwd()
# p = Path(path)

# concatfile = '/home/meiker/phylo_tree/iqtree/concat_alignments'
# #concat_file = '/home/meike/tests/Files/concat_alignments'
# #testfile = '/home/meike/tests/Files/prottest_testset1.fa'

# seqs = {}
# with open(concatfile) as f:
#     name = ''
#     for line in f:
#         line = line.strip()
#         if line.startswith('>'):
#             name = line[1:]
#             seqs[name] = ''
#         else:
#             seqs[name] += line
            


# # generate a map of unique sequence to identical sequences
# similar = []

# names = []
# sequences = []
# for name, seq in seqs.items():
#     names.append(name)
#     sequences.append(seq)

# #ratios = []
    
# #iterate over the sequences and compare each sequence with each other
# for a, b in itertools.combinations(sequences, 2):
#     Ratio = fuzz.ratio(a,b)
#     #print(count, ' Ratio determined')
    
#     #if the sequences are at least 95% identical add the ids (names) to the list 'similar'
#     if Ratio >= 95:
#         #ratios.append(Ratio)
#         a_ind = sequences.index(a)
#         b_ind = sequences.index(b)
#         similar.append((names[a_ind], names[b_ind]))

# #make dict to store ids (first appeared) as keys and all ids with similar seqs as values
# ids2dupids = {}

# left_over_ids = []
# for i in similar:
#     #check if id is already a key or values (list of duplicate ids)
#     if i[0] not in ids2dupids.keys() and i[0] not in [x for v in ids2dupids.values() for x in v]:
#         ids2dupids[i[0]] = [i[1]]
#         left_over_ids.append(i[0])
#     else:
#         #loop through dict and look at keys and duplicates
#         for key, duplicate in ids2dupids.items():
#             #if the first item of the similar sequences is a key, than add the second to the values as duplicate
#             if i[0] == key:
#                 ids2dupids[key].append(i[1])
# print("Number of left over IDs: ", len(left_over_ids),"\nleft over IDs: ", left_over_ids)
                

# unique_seq_ids = ids2dupids.keys()
              
# with open('/home/meiker/phylo_tree/iqtree/reduced_alignments/reduced_concat_alignments.fa', 'w') as f:
#     for k in unique_seq_ids:
#         samp = (k, seqs[k])
#         f.write(">{}\n{}\n".format(*samp))
            

            
