#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 16:31:04 2019

@author: meike
"""
'''Test to summarize and Visualize the metadata --> gneome lengths'''

import matplotlib.pyplot as plt
import os
from pathlib import Path

path = os.getcwd()
p = Path(path)


def size_grouping(sizes):  
    """
    list of all sizes. Groups intp ranges of 5000,10000,15000 and 20000+
    """
    lengths = [i for i in range(round(min(sizes), -5), round(max(sizes), -5), 100000)]
    
    groups = {k: [] for k in lengths}
    
    for i in sizes:
        counted = False
        for group in lengths:
            if i in range(round(group - 50000), round(group + 50000)) and counted == False:
                groups[group].append(i)
                counted = True
    
    return groups

def size_count(sizes):  
    """
    list of all sizes. Groups int ranges of 1,700,000 in steps of 100,000 up t0 2,800,000.
    counts the amount of species within a range of plus/min 50,000 of beforementioned ranges.
    """
    lengths = [i for i in range(round(min(sizes), -5), round(max(sizes), -5), 100000)]
    
    groups = {k: 0 for k in lengths}
    
    for i in sizes:
        counted = False
        for group in lengths:
            if i in range(round(group - 50000), round(group + 50000)) and counted == False:
                groups[group] += 1
                counted = True
    
    return groups

def genome_size_plot(inputfile, saving_dir, genus):
    '''
    Needs inputfile containing column with genomes lengths. Groups genomes sizes in 
    ranges from min to max with steps of 100,000 and save a tiff of a bar plot in given dir.
    '''
    genus = genus.capitalize()
    sizes = []
    with open (inputfile) as f:
        headers = f.readline().strip().split('\t')
        headers_inds = {name: i for i, name in enumerate(headers)}
        
        for line in f:
            a = line.strip().split('\t')
            for i, item in enumerate(a):
                if headers_inds['genome_length'] == i and item != '':
                    sizes.append(int(item))
        counts = size_count(sizes)
        
        labels =[]
        for k in counts:
            labels.append(k/1000000)
            
        plt.bar(range(len(counts)), list(counts.values()), align='center')
        plt.xticks(range(len(counts)), labels, rotation=90)
        plt.title("Genome sizes " + genus, fontsize = 14)
        plt.xlabel("Genome length\n(in Mb)", fontsize = 11)
        plt.ylabel("Number of species counted", fontsize = 11)
        plt.tight_layout()
        plt.savefig(saving_dir)

lactofile = os.path.join(p.parents[0], 'files', '06012020_lactococcus_database.tsv')
lacto_save = os.path.join(p.parents[0], 'figures', '20200107_lactococcus_genome_lengths.png')

genome_size_plot(lactofile, lacto_save, 'lactococcus')

florifile = os.path.join(p.parents[0], 'files', '06012020_floricoccus_database.tsv')
flori_save = os.path.join(p.parents[0], 'figures', '20200107_floricoccus_genome_lengths.png')
genome_size_plot(florifile, flori_save, 'floricoccus')


streptofile = os.path.join(p.parents[0], 'files', '06012020_streptococcus_database.tsv')
strepto_save = os.path.join(p.parents[0], 'figures', '20200107_streptococcus_genome_lengths.png')
genome_size_plot(streptofile, strepto_save, 'streptococcus')



