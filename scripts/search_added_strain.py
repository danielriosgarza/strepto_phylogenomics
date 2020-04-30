#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 10:54:46 2020

@author: meike
"""
a = '/home/meiker/phylo_tree/roary/29042020_reduced_concat_alignments.fa'
b = '/home/meiker/phylo_tree/iqtree/reduced_alignments/09032020_reduced_concat_alignments.fa'

l = []
with open(a) as f:
    with open(b) as f2:
        for line in f:
            if line.startswith('>'):
                l.append(line)
        for line in f2:
            if line.startswith('>'):
                if line not in l:
                    print(line)
