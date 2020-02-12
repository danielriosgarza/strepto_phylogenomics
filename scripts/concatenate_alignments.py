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
    


sequences={}

path='/home/meiker/phylo_tree/msa_trimmed/'

files = os.listdir(path)

for fl in files[0:2]:
    with open(path+fl) as f:
        f.readline()
        f.readline()
        for line in f:
            if line!='\n':
                a= line.strip()
                seq = a[24::]
                id_ = a[0:24].split('_')
                id_ = id_[0][0:5] + id_[1][0:5]
                
                if id_ not in sequences:
                    sequences[id_]=''
                sequences[id_]+=seq

write_phylip(sequences, '/home/meiker/phylo_tree/concat_alignment.phy')


