#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 10:04:37 2020

@author: meike
"""


import os


def write_phylip(seqD, filePath):
    seqIds = list(seqD.keys())
    #get len of first seqence
    l1 = len(seqD[seqIds[0]])
    
    for i in seqIds:
        if len(seqD[i])!=l1:
            print('Error: Alignment of different length')
            return None
    
    with open(filePath, 'w') as f:
        f.write(str(len(seqIds))+' ' + str(l1) +'\n')
        
        for seq in seqIds:
            f.write(seq + seqD[seq]+'\n')

            
all_alignments ={}

path='/home/meike/strepto_phylogenomics/files/random_files/test_phylip/'

t_file ='/home/meike/strepto_phylogenomics/files/random_files/test_for_phylipformat'

files = os.listdir(path)

ids = ['streptococcus_0000'+str(i) for i in range(1,6)]


for fl in files:
    sequences ={}
    aligned_ids = []
    with open(path+fl) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                seq=''
                a = line.split('_')
                phy_id = a[0][1:6] + a[1]
                db_id = line[1::]
                if db_id not in aligned_ids:
                    aligned_ids.append(db_id)
                if phy_id not in sequences:
                    sequences[phy_id]=''
            else:
                seq+=line
            sequences[phy_id]=seq
    
    for id_ in ids:
        if id_ not in aligned_ids:
            a = id_.split('_')
            phy_id = a[0][0:5] + a[1]
            sequences[phy_id] =''
            
    longest = len(max(sequences.values()))
    
    for k, v in sequences.items():
        if len(v) < longest:
            difference = longest-len(v)
            gap = '-'*difference
            sequences[k] = v+gap
    all_alignments[fl] = sequences        

with open('/home/meike/strepto_phylogenomics/files/random_files/concatenated_test', 'w') as f:
    for gene in all_alignments:
        seqD = all_alignments[gene] #get dict with alignments per gene
        number_rows = str(len(seqD.keys())) #determine number of ids (first line phylip format)
        seq_length = str(len(list(seqD.values())[0])) #determine alignment length (first line phylip format)
        f.write(number_rows+ ' ' + seq_length+ '\n')
        for id_, seq in seqD.items():
            f.write(id_+seq+'\n')
            



#write_phylip(sequences,'/home/meike/strepto_phylogenomics/files/random_files/phylip_format.phy')


#%%
        
files = os.listdir(path)

path='/home/meiker/phylo_tree/msa_test/msa_trimmed/'

#make dict containing ids for iqtree (10 chars) as keys and alignment as sequences
sequences={}


for fl in files:
    with open(path+fl) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                seq=''
                a = line.split('_')
                phy_id = a[0][1:6] + a[1]
                db_id = line[1::]
                if db_id not in aligned_ids:
                    aligned_ids.append(db_id)
                if phy_id not in sequences:
                    sequences[phy_id]=''
            else:
                seq+=line
            sequences[phy_id]=seq
        aligned_ids = []
        for id_ in ids:
            if id_ not in aligned_ids:
                a = id_.split('_')
                phy_id = a[0][0:5] + a[1]
                sequences[phy_id] =''
        
longest = len(max(sequences.values()))

for k, v in sequences.items():
    if len(v) < longest:
        difference = longest-len(v)
        gap = '-'*difference
        sequences[k] = v+gap

write_phylip(sequences, '/home/meiker/phylo_tree/concat_alignment.phy')



